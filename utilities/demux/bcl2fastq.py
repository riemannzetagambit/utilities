#!/usr/bin/env python

import argparse
import glob
import os
import re
import subprocess
import sys

import pandas as pd

from utilities.log_util import get_logger, log_command


BCL2FASTQ = 'bcl2fastq'

S3_RETRY = 5
S3_LOG_DIR = 's3://trace-genomics/TraceGenomics/logs'
ROOT_DIR_PATH = '/tmp'


def get_default_requirements():
    return argparse.Namespace(vcpus=64, memory=256000, storage=2000,
                              queue='aegea_batch_demux',
                              ecr_image='demuxer',
                              ulimits=['nofile:1000000'])


def get_parser():
    parser = argparse.ArgumentParser(
            prog='bcl2fastq.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--exp_id',
                        help='Name of the prefix on S3 in --s3_input_dir that contains the run you want to process',
                        required=True)

    parser.add_argument('--s3_input_dir',
                        default='s3://trace-genomics/TraceGenomics',
                        help='S3 path for [exp_id] folder of BCL files')
    parser.add_argument('--s3_output_dir',
                        default='s3://trace-genomics/TraceGenomics',
                        help='S3 path to put fastq files')
    parser.add_argument('--s3_report_dir',
                        default='s3://trace-genomics/TraceGenomics/reports',
                        help='S3 path to put the bcl2fastq report')
    parser.add_argument('--s3_sample_sheet_dir',
                        default='s3://trace-genomics/TraceGenomics/sample-sheets',
                        help='S3 path to look for the sample sheet')

    parser.add_argument('--group_by_sample', action='store_true',
                        help='Group the fastq files into folders based on sample name')
    parser.add_argument('--skip_undetermined', action='store_true',
                        help="Don't upload the Undetermined files (can save time)")
    parser.add_argument('--no_s3_download', action='store_true',
                        help="Do not download bcl files from S3 (useful if testing or already have locally "
                             "demultiplexed files in the ROOT_DIR_PATH location. Currently a work in progress.")
    parser.add_argument('--no_s3_upload', action='store_true',
                        help="Do not upload demultiplexed fastq files to S3 (useful if testing and don't want to "
                             "check demultiplexing locally without writing to S3. Currently a work in progress.")

    parser.add_argument('--sample_sheet_name', default=None,
                        help='Defaults to [exp_id].csv')
    parser.add_argument('--force-glacier', action='store_true',
                        help='Force a transfer from Glacier storage')
    # TODO(dstone): add an option to delete the original un-demultiplexed from S3 afterward

    parser.add_argument('--bcl2fastq_options',
                        default=['--no-lane-splitting'],
                        nargs=argparse.REMAINDER,
                        help='Options to pass to bcl2fastq')

    return parser


def _check_for_run_information(sample_name):
    """
    Helper function to try and find run ID information of the form RunXX_YY
    """
    m = re.match('^Run\d+_\d+$', sample_name)
    if m is not None:
        return True
    else:
        return False


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get('AWS_BATCH_JOB_ID'):
        root_dir = os.path.join(ROOT_DIR_PATH, os.environ['AWS_BATCH_JOB_ID'])
    else:
        root_dir = ROOT_DIR_PATH

    if args.sample_sheet_name is None:
        args.sample_sheet_name = '{}.csv'.format(args.exp_id)

    # local directories
    result_path = os.path.join(root_dir, 'data', 'hca', args.exp_id)
    bcl_path = os.path.join(result_path, 'bcl')
    output_path = os.path.join(result_path, 'fastqs')

    if not args.no_s3_download:
        # only make dirs if they don't exist yet
        if not os.path.isdir(result_path):
            os.makedirs(result_path)
        if not os.path.isdir(bcl_path):
            os.mkdir(bcl_path)

        # download sample sheet
        command = ['aws', 's3', 'cp', '--quiet',
                   os.path.join(args.s3_sample_sheet_dir, args.sample_sheet_name),
                   result_path]
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying s3 copy")
        else:
            raise RuntimeError("couldn't download sample sheet {}".format(
                    os.path.join(args.s3_sample_sheet_dir, args.sample_sheet_name))
            )

        # do a check on the sample inputs to make sure we can get run IDs from all of them
        # change this if the Illumina sample sheet output ever changes; otherwise this line has the headers
        _SAMPLE_SHEET_STARTING_LINE = 21
        df_csv = pd.read_csv(os.path.join(result_path, args.sample_sheet_name), header=_SAMPLE_SHEET_STARTING_LINE)
        samples_not_matching_run_ids = [sample_name for sample_name in df_csv['Sample_ID'] if not _check_for_run_information(sample_name)]
        if len(samples_not_matching_run_ids) > 0:
            raise ValueError('Found sample names that I could not extract run ID values (of the form RunXX_YY) from: '
                             '{}'.format(samples_not_matching_run_ids))

        # download the bcl files
        command = ['aws', 's3', 'sync', '--quiet',
                   '--force-glacier-transfer' if args.force_glacier else '',
                   os.path.join(args.s3_input_dir, args.exp_id), bcl_path]
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying s3 sync bcl")
        else:
            raise RuntimeError("couldn't sync {}".format(
                    os.path.join(args.s3_input_dir, args.exp_id))
            )


    # this is actually awful because the process forks and you have to go kill it yourself
    command = ('while true;'
               ' do memusage=`cat /sys/fs/cgroup/memory/memory.usage_in_bytes`;'
               ' memgb=`echo "${memusage}/(1000000000)" | bc -l | xargs -I {} printf "%.2f\n" {}`;'
               ' echo "memory usage: ${memgb}GB";'
               ' echo "disk usage: " `df -h | grep -e "/$" | awk \'{print $(NF-4)" "$(NF-3)" "$(NF-2)" "$(NF-1)" "$NF}\''
               ' sleep 90;'
               ' done')
    p = subprocess.Popen([command], shell=True)

    # Run bcl2 fastq
    command = [BCL2FASTQ, ' '.join(args.bcl2fastq_options),
               '--sample-sheet', os.path.join(result_path,
                                              args.sample_sheet_name),
               '-R', bcl_path, '-o', output_path]
    log_command(logger, command, shell=True)

    # fix directory structure of the files *before* sync!
    fastqgz_files = glob.glob(os.path.join(output_path, '*fastq.gz'))
    logger.debug('all fastq.gz files\n{}\n\n'.format('\n'.join(fastqgz_files)))

    # TODO(dstone): organize the run based on the TraceGenomics/RunXX/RunXX_YY/*.fastq.gz and do our usual rearrangement
    for fastq_file in fastqgz_files:
        if (args.skip_undetermined
            and os.path.basename(fastq_file).startswith('Undetermined')):
            logger.info("removing {}".format(os.path.basename(fastq_file)))
            os.remove(fastq_file)
        elif args.group_by_sample:
            # exclude the sample number (_S[numbers])
            m = re.match("(.+)(_S\d+_R[12]_001.fastq.gz)",
                         os.path.basename(fastq_file))
            if m:
                sample = m.group(1) # should be of the form RunX_Y
                if not re.match('^Run\d+_\d+$', sample):
                    # shouldn't actually be able to get here, because there is a check above at the sample sheet level,
                    # but just in case
                    raise ValueError('Was expecting to find a sample name of the form RunXX_YY, could not find in {} sample name!'.format(sample))
                run = sample.split('_')[0]
                grouped_sample_path = os.path.join(output_path, run, sample) # organizes as RunX/RunX_Y/[sample stuff]
                if not os.path.exists(grouped_sample_path):
                    logger.debug("creating {}".format(grouped_sample_path))
                    os.makedirs(grouped_sample_path)
                logger.debug("moving {}".format(fastq_file))
                os.rename(fastq_file, os.path.join(grouped_sample_path, os.path.basename(fastq_file)))
            else:
                logger.warning("Warning: regex didn't match {}".format(fastq_file))

    sys.stdout.flush()

    if not args.no_s3_upload:
        # upload fastq files to destination folder
        command = ['aws', 's3', 'sync', '--quiet', output_path,
                   args.s3_output_dir,
                   # this doesn't fit our output structure
                   #os.path.join(args.s3_output_dir, args.exp_id, 'rawdata'),
                   '--exclude', '"*"', '--include', '"*fastq.gz"']
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying sync fastq")
        else:
            raise RuntimeError("couldn't sync fastqs")


        # check fastq upload
        command = ['aws', 's3', 'ls', '--recursive',
                   args.s3_output_dir]
                   #os.path.join(args.s3_output_dir, args.exp_id, 'rawdata')]
        log_command(logger, command, shell=True)


        # Move reports data back to S3
        reports_path = subprocess.check_output(
                "ls -d {}".format(os.path.join(output_path, 'Reports', 'html', '*',
                                               'all', 'all', 'all')),
                shell=True).rstrip()
        command = ['aws', 's3', 'cp', '--quiet', reports_path,
                   os.path.join(args.s3_report_dir, args.exp_id),
                   '--recursive']
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying cp reports")
        else:
            raise RuntimeError("couldn't cp reports")

    p.kill()


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        if log_file:
            log_cmd = 'aws s3 cp --quiet {} {}'.format(log_file, S3_LOG_DIR)
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
