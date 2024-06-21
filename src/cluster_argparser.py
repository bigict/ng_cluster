import argparse
import textwrap
import os
from split import split_parameters_file
from aggregation import aggregate_result_file

script_dir = os.path.dirname(os.path.realpath(__file__))

def get_usage_info():
    # read in the example_commands.txt file and return it as a string
    f = open(os.path.join(script_dir, 'example_commands.txt'), 'r')
    lines = f.readlines()
    return "".join(lines)

class ClusterArgumentParser:
    parser = argparse.ArgumentParser(
        usage=textwrap.dedent(get_usage_info()),
        description=textwrap.dedent('''
            This tool groups epitopes into clusters based on sequence identity. 
            A cluster is defined as a group of sequences which have a sequence similarity 
            greater than the minimum sequence identity threshold specified.
            '''),
        formatter_class=argparse.RawTextHelpFormatter)

    def __init__(self):
        pass

    def parse_args(self):
        
        #TODO: somehow argparse is calling json-input optional parameters in the help message, even though it is marked as required
        self.parser.add_argument("--json-input", "-j",
                                 dest="json_filename",
                                 help="JSON file containing input parameters.",
                                 metavar="JSON_FILE")
        self.parser.add_argument("--output-prefix", "-o",
                                 dest="output_prefix",
                                 help="prediction result output prefix.",
                                 metavar="OUTPUT_PREFIX")
        self.parser.add_argument("--output-format", "-f",
                                 dest="output_format",
                                 default="tsv",
                                 help="prediction result output format (Default=tsv)",
                                 metavar="OUTPUT_FORMAT")        
        
        # TODO: consider hiding most of these
        # Optional Arguments (Flags)
        self.parser.add_argument("--split",
                                 action="store_true",
                                 dest="split_parameters_flag",
                                 default=False,
                                #help="flag to indicate the action we want to take with the standalone: split parameters into JSON files")
                                help=argparse.SUPPRESS)
        self.parser.add_argument("--split-dir",
                                 dest="split_parameters_dir",
                                 default='',
                                 #help="the directory for the JSON files that input parameters are split into")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--split-inputs-dir",
                                 dest="split_inputs_dir",
                                 default=None,
                                 #help="the directory for the sequence and peptide files that input sequences split into")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--aggregate",
                                 action="store_true",
                                 dest="aggregate_parameters_flag",
                                 default=False,
                                 #help="flag to indicate the action to aggregate the results")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--job-desc-file",
                                 dest="job_desc_file",
                                 default='',
                                 #help="the file path for the job description")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--aggregate-input-dir",
                                 dest="aggregate_input_dir",
                                 default='',
                                 #help="the directory for the JSON files which have input parameters")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--aggregate-result-dir",
                                 dest="aggregate_result_dir",
                                 default='',
                                 #help="the directory for the JSON files contains results need to be aggregated as well as the place we place the final result file")
                                 help=argparse.SUPPRESS)
        self.parser.add_argument("--assume-valid",
                                 action="store_true",
                                 dest="assume_valid_flag",
                                 default=False,
                                 #help="flag to indicate skipping validation")
                                 help=argparse.SUPPRESS)
        
        return self.parser.parse_args(), self.parser


    def aggregate_parameters(self, args) :
        if getattr(args, 'aggregate_parameters_flag') :
            job_desc_file = getattr(args, 'job_desc_file')
            aggregate_input_dir = getattr(args, 'aggregate_input_dir')
            aggregate_result_dir = getattr(args, 'aggregate_result_dir')
            aggregate_result_file(job_desc_file, aggregate_input_dir, aggregate_result_dir)
            exit(0)

    def split_parameters(self, args) :
        if getattr(args, 'split_parameters_flag') :
            json_filename = getattr(args, 'json_filename')
            split_parameters_dir = getattr(args, 'split_parameters_dir')
            split_inputs_dir = getattr(args, 'split_inputs_dir')
            assume_valid_flag = getattr(args, 'assume_valid_flag')
            split_parameters_file(json_filename, split_parameters_dir, split_inputs_dir, assume_valid=assume_valid_flag)
            exit(0)
