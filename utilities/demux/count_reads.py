import argparse

import boto3


def arguments():
    parser = argparse.ArgumentParser(description='Process hits in our DB to indicators for the application database.')
    parser.add_argument('-run',
                        type=int,
                        help='Run number to process. Example: for Run101, use "101".',
			required=True)
    parser.add_argument('-bucket',
                        type=str,
                        help='Name of S3 bucket to read from.',
                        default='trace-genomics')
    parser.add_argument('-logname',
                        type=str,
                        help='Name of logfile to use.',
                        default=None)
    args = parser.parse_args()
    return args, parser


def _get_total_reads_for_run(run, s3_bucket='trace-genomics'):
    # should probably abstract this at some point and do something like du.get_database_connection() with boto clients to enforce singleton patterns
    s3 = boto3.client('s3')

    reads_dict = {}
    # don't duplicate read counts, only count reads from forward fastq file
    _FASTQ_FORWARD_REGEX = re.compile('^(?P<run_id>Run\d+_\d+)_S\d+_R1_001.fastq.gz$')

    s3lo = s3.list_objects(Bucket=s3_bucket, Prefix='sequencing_data/TraceGenomics/Run{}/'.format(run))
    for o in s3lo['Contents']:
        m = _FASTQ_FORWARD_REGEX.match(o['Key'].split('/')[-1])
        if m:
            run_id = m.group('run_id')
            r = s3.select_object_content(
                Bucket=s3_bucket,
                Key=o['Key'],
                ExpressionType='SQL',
                # this line counts the number of reads in the file
                Expression="select count(*) from s3object s",
                # because we delimit on @, which is unique for each 4 lines/read, this line insures that the count above is the number of reads (not number of lines)
                InputSerialization = {'CSV': {'RecordDelimiter': '@', 'FieldDelimiter': '\n'}, 'CompressionType': 'GZIP',},
                OutputSerialization = {'CSV': {}},)

            # assume that the 'Records' has a subdict with key 'Payloads';
            # the number returned gives you the total number of lines/4 in the fastq.gz file of interest, which is equivalent to the total number of reads
            # the 'RecordDelimiter' kwarg in InputSerialization above gives us each set of 4 lines as one record and does the lines -> reads work for us
            total_reads = [float(event['Records']['Payload'].strip()) for event in filter(lambda x: x.get('Records') is not None, r['Payload'])]

            # furthermore assume that there is only one event that has the above information
            if len(total_reads) != 1:
                raise ValueError('Expected only one measurement of total lines in fastq file, instead got {} possible values return from boto EventStream'.format(len(total_lines)))
            total_reads = total_reads[0]
            reads_dict[run_id] = total_reads
    return reads_dict


def _output_reads_for_run_to_csv(run, output_file=None, **kwargs):
    reads_dict = _get_total_reads_for_run(run, **kwargs)
    # use pandas as a convenience to write to csv for us
    # clean up and sort the run IDs in the 'natural' way with sort_for_humans
    df_reads = pd.DataFrame.from_dict(reads_dict, orient='index')
    df_reads = df_reads.loc[sort_for_humans(df_reads.index)]

    if output_file is None:
        output_file = './run{}_reads.csv'.format(run)
    df_reads.to_csv(output_file, header=None)


def main(args):
    _output_reads_for_run_to_csv(args.run)


if __name__ == '__main__':
    args, _  = arguments()
    if args.logname is None:
        logname = __file__.replace('.py', '.log')
        logname = os.path.abspath(logname)
    else:
        logname = args.logname
    logname = os.path.abspath(logname)

    # goes to the console and will ignore further basicConfig calls. Remove the handler if there is one.
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    # this appears to do nothing in python 2
    logging.basicConfig(filename=logname,
                        filemode='w',
                        format='%(asctime)s %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging_level)
    baselog = logging.getLogger(__name__)
    baselog.setLevel(logging_level)

    main(args)
