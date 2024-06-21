
#!/usr/bin/env python
import os
import sys
import re
import tempfile
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', datefmt='%Y-%m-%d:%H:%M:%S',)

from epitope_cluster_analysis.analysis import process

def main():
    import select

    try:
        usage = ""

        #TODO: Add option to run in 'internal mode' that will allow
        #      skipping certain validations, e.g.:
        #      checking if peptides are valid, method/length, etc.

        parser = OptionParser(usage=usage, version="%prog {}".format(self.version), add_help_option=False)
        parser.add_option("-v", "--versions",
                          action="store_true",
                          dest="version_flag",
                          default=False,
                          help="print specific methods and their versions.")
        parser.add_option("-h", "--help",
                          action="store_true",
                          dest="help",
                          default=False,
                          help="print available commands.")

        parser.add_option("-f", dest="seqfile_name",
                          help="FILE containing a peptide sequences in fasta format.", metavar="FILE")

        parser.add_option("-t", dest="threshold",
                          help="threshold.", metavar="THRESHOLD")

        parser.add_option("--method", dest="method",
                          help="prediction method.", metavar="METHOD")

        parser.add_option("--min", dest="minimum_length", default='1',
                          help="minimum_length", metavar="minimum_length")
        parser.add_option("--max", dest="maximum_length", default='100',
                          help="maximum_length", metavar="MAXIMUM_LENGTH")

        parser.add_option("--output", "-o", dest="output",
                          help="prediction result putput path.", metavar="OUTPUT")

        (options, args) = parser.parse_args()

        if options.help:
            commandline_help()
            exit(0)

        if not sys.stdin.isatty():
            stdin = sys.stdin.readline().strip()
            args.append(stdin)

        args = list(filter(None, args))
        commandline_input_prediction(options, args)


    except UnexpectedInputError as e:
        print( str(e))



def print_result(result):
    print('Cluster.Sub-Cluster Number\tPeptide Number\tAlignment\tPosition\tDescription\tPeptide\tCluster Consensus')
    table_rows = result
    for row in table_rows:
        print( '\t'.join(map(str, row)))

def commandline_help():
    return

def commandline_input_prediction(options, args):
    """ This version takes a file containing an peptide sequences as input."""

    # 1. read input params
    seqfile_name = options.seqfile_name
    method = options.method
    minimum_length = options.minimum_length
    maximum_length = options.maximum_length
    threshold = options.threshold

    # 2 validation
    errors =  input_validation(seqfile_name, threshold, minimum_length, maximum_length, method)
    if errors:
        print(errors)

    # 3. predict
    result = process(seqfile_name, threshold, minimum_length, maximum_length, method)

    # 4. output
    print_result(result)



def input_validation(seqfile_name, threshold, minimum_length, maximum_length, method):
    return


if __name__ == '__main__':
    main()
