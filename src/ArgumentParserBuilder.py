import sys
import warnings
import argparse
import textwrap
# import ArgumentValidators as arg_validator
from pathlib import Path
from abc import ABC


class NGArgumentParserBuilder(ABC):
    def __init__(self):
        ''' NOTE
        
        Available options can be found in the argparse document:
        https://docs.python.org/3/library/argparse.html#action

        TIPS:
        nargs - '?' : 'default' value is used when no command-line arg is present.
        'const' value is used when the flag is listed, but no argument is provided.
        ex) parser.add_argument('--foo', nargs='?', const='c', default='d')
            >> parser.parse_args(['XX', '--foo'])
            >> Namespace(bar='XX', foo='c')
            >> parser.parse_args([])
            >> Namespace(bar='d', foo='d')
        '''
        self.parser = argparse.ArgumentParser()


    ###########################################################
    # Setting Program Descriptions
    ###########################################################
    @property
    def name(self) -> str:
        return self.parser.name
    
    @name.setter
    def name(self, name: str='PROG') -> None:
        self.parser.prog = name 

    @property
    def description(self) -> str:
        return self.parser.description

    @description.setter
    def description(self, text: str) -> None:
        self.parser.description = text
        self.parser.formatter_class = argparse.RawDescriptionHelpFormatter

    @property
    def epilog(self) -> str:
        return self.parser.epilog
    
    @epilog.setter
    def epilog(self, epilog: str) -> None:
        self.parser.epilog = epilog

    @property
    def version(self) -> str:
        return self._version
    
    @version.setter
    def version(self, version: str) -> None:
        self.parser.add_argument(
            '-v',
            '--version',
            action='version',
            version=f'%(prog)s {version}'
        )

        self._version = version

    def parse_arguments(self) -> argparse.Namespace:
        return self.parser.parse_args()

    ###########################################################
    # Predefined Positional Arguments
    ###########################################################
    def sequences(self, 
                  action: object=None, 
                  required: bool=False) -> None:
        '''
        Positional argument that will take space separated sequences
        as list of strings.
        '''
        narg_value = '*'
        custom_action = 'store'

        # If the tool absolutely requires sequences and would like
        # to throw error message if none are provided, use '+'.
        if required: narg_value = '+'

        if action:
            custom_action = action

        self.parser.add_argument(
            'sequences',
            type=str,
            nargs=narg_value,
            action=custom_action,
            help='Space separated peptide sequences.')

    ###########################################################
    # Predefined Optional Arguments
    ###########################################################
    def alleles(self, 
                tool: str='mhci',
                action: object=None) -> None:
        ''' Rebranding from allele-list
        
        Originally 'allele-list' is used, which took list of alleles.
        Simply renaming it to 'alleles' as the name implies multiple alleles.
        '''
        
        # Add 'tool' argument quietly to the parser.
        # Suppressing help, so that this is option is hidden away
        # from users.
        # NOTE:
        # This will variable will be used in ValidateAlleles
        self.parser.add_argument(
            'tool',
            nargs='?',
            type=str,
            default=tool,
            const=tool,
            help=argparse.SUPPRESS
        )

        self.parser.add_argument(
            '-a', 
            '--alleles', 
            type=str,
            # Gathers all arguments as a list
            nargs='*', 
            # action=arg_validator.ValidateAlleles,
            help="List of alleles."
        )

    # TODO
    def lengths(self) -> None:
        '''peptide length'''
        pass

    
    # NOTE: What to do with version of the method?
    def method(self) -> None:
        pass

    # TODO
    def threshold(self) -> None:
        pass

        
    def database_path(self, 
                      single_dash_option: str, 
                      double_dash_option: str, 
                      help: str, 
                      action: object=None, 
                      required: bool=False, 
                      dest: str='db_path') -> None:
        help_txt = '''
            Sets path to the database.
            '''
        
        # 'action' default -- according to python doc
        custom_action = 'store'
        
        if help: help_txt = help

        if action: 
            custom_action = action
        else:
            warnings.warn('''\
                \nThe \'database_path\' option has not been validated.\
                \nPlease implement custom validation.''', 
                UserWarning
            )

        # Allow users to set path to the database
        self.parser.add_argument(
            f'-{single_dash_option}',
            f'--{double_dash_option}',  
            dest = dest, 
            required = required,
            nargs = '?',
            type = Path,
            const = None,
            action = custom_action,
            help = textwrap.dedent(help_txt)
        )
        


    # Types of files :
    #     FASTA
    #     CSV
    #     TSV
    #     TXT - list of sequences
    def input_file(self, 
                   single_dash_option: str='i',
                   double_dash_option: str='input-file',
                   dest: str='',
                   action: object=None,
                   required: bool=False,
                   help: str='') -> None:
        custom_action = 'store'
        help_text = 'a TXT/CSV/FASTA file containing protein sequences.'
        
        if help: help_text = help
        
        # Default behavior of 'dest' if not given.
        if not dest:
            dest = double_dash_option.replace('-', '_')

        if action:
            custom_action = action


        self.parser.add_argument(
            f'-{single_dash_option}',
            f'--{double_dash_option}',
            dest = dest, 
            required = required,
            nargs = '?', 
            type = argparse.FileType('r'),
            # NOTE: File type (csv, fasta, tsv) validation can be set here.
            action=custom_action,
            # Do not add this attribute to Namespace() if not provided
            default = argparse.SUPPRESS, 
            help = help_text
        )


    # Types of output file:
        # JSON
        # TSV
        # CSV
    def output_file(self, 
                    single_dash_option: str='o',
                    double_dash_option: str='output-file',
                    dest: str='',
                    action: object=None,
                    required: bool=False,
                    help: str='') -> None:
        custom_action = 'store'
        help_text = 'Output file containing the prediction result table.'

        if help: help_text = help

        if not dest:
            dest = double_dash_option.replace('-', '_')
        
        if action:
            custom_action = action

        self.parser.add_argument(
            f'-{single_dash_option}',
            f'--{double_dash_option}',
            dest = dest, 
            required = required,
            nargs = '?', 
            # Type is set from 'argparse.FileType('w') to 'str' as it needs to 
            # combine wih 'output_format' to get the full file name. It will 
            # then be open with 'with' keyword.
            type = str,
            # NOTE: File type (json, csv, tsv) validation can be set here.
            action=custom_action,
            # Do not add this attribute to Namespace() if not provided
            default = argparse.SUPPRESS, 
            help = help_text
        )

    def output_format(
            self, 
            single_dash_option: str='f',
            double_dash_option: str='output-format',
            dest: str='',
            choices: list=None,
            default: list=None,
            action: object=None,
            required: bool=False,
            help: str='') -> None:
        
        custom_action = 'store'
        help_text = 'Output file containing the prediction result table.'

        if help: help_text = help

        if not dest:
            dest = double_dash_option.replace('-', '_')
        
        if action:
            custom_action = action

        
        self.parser.add_argument(
            f'-{single_dash_option}',
            f'--{double_dash_option}',
            dest = dest, 
            required = required,
            nargs = 1,
            type = str,
            choices = choices,
            action = custom_action,
            default = default,
            help = help_text
        )

    
    def data_dir(self) -> None:
        self.parser.add_argument(
            "--data-dir",
            type=Path,
            default=Path(__file__).absolute().parent / "data",
            help="Path to the data directory",
        )


    ###########################################################
    # Split and aggregate flags.
    ###########################################################
    # NOTE:
    # These flags will be less flexible compared to other flags.
    def split(self) -> None:
        self.parser.add_argument(
            "--split",
            action="store_true",
            dest="split_parameters_flag",
            default=False,
            help="flag to indicate the action we want to take with the standalone: split parameters into JSON files"
        )

    def split_dir(self) -> None:
        self.parser.add_argument(
            "--split-dir",
            dest="split_parameters_dir",
            default='',
            help="the diretory for the JSON files that input parameters splitted into"
        )

    def split_input_dir(self) -> None:
        self.parser.add_argument(
            "--split-inputs-dir",
            dest="split_inputs_dir",
            default=None,
            help="the diretory for the sequence and peptide files that input sequences splitted into"
        )

    def aggregate(self) -> None:
        self.parser.add_argument(
            "--aggregate",
            action="store_true",
            dest="aggregate_parameters_flag",
            default=False,
            help="flag to indicate the action to aggregate the results"
        )

    def aggregate_input_dir(self) -> None:
        self.parser.add_argument(
            "--aggregate-input-dir",
            dest="aggregate_input_dir",
            default='',
            help="the diretory for the JSON files which have input parameters"
        )

    def aggregate_result_dir(self) -> None:
        self.parser.add_argument(
            "--aggregate-result-dir",
            dest="aggregate_result_dir",
            default='',
            help="the diretory for the JSON files contains results need to be aggregated as well as the place we place the final result file"
        )
    
    def aggregate_output_format(self,
                                file_type: str='tsv') -> None:
        self.parser.add_argument(
            "-af",
            "--aggregate-output-format",  
            dest="aggregate_output_format", 
            default=file_type,
            help="prediction result output format.", 
            metavar="OUTPUT_FORMAT")


    def job_description_file(self) -> None:
        self.parser.add_argument(
            "--job-desc-file",
            dest="job_desc_file",
            default='',
            help="the file path for the job description"
        )
    
    def assume_valid(self) -> None:
        self.parser.add_argument(
            "--assume-valid",
            action="store_true",
            dest="assume_valid_flag",
            default=False,
            help="flag to indicate skiping validation"
        )

    def json_filename(self) -> None:
        self.parser.add_argument(
            "-j",
            dest="json_filename",
            required=False,
            help="JSON file containing all parameters.",
            metavar="JSON_FILE")

