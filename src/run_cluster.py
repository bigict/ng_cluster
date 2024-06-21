import os
import sys
import json
import tempfile
import logging
from urllib.request import urlopen
from shutil import copyfileobj
from cluster_argparser import ClusterArgumentParser
from validations import cluster_validate
from predictions import cluster_analysis

import warnings

# Suppress all warnings
warnings.filterwarnings("ignore")

# logging config
logging.basicConfig(level=logging.INFO, format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', datefmt='%Y-%m-%d:%H:%M:%S',)

# adding all methods to the python path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir )
methods_dir = os.path.join(script_dir, '../method')

for method_dir_name in os.listdir(methods_dir):
    method_base_dir = os.path.join(methods_dir, method_dir_name)
    if os.path.isdir(method_base_dir):
        sys.path.append(method_base_dir)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def truncate_file(output_path):
    with open(output_path, 'a+') as tsv_file:
        tsv_file.truncate(0)

def print_result(result):
    table_rows = result
    for row in table_rows:
        print( '\t'.join(map(str, row)))

def save_tsv(result, output_path):
    table_rows = result
    with open(output_path, 'a') as tsv_file:
        tsv_file.write( '\n'.join(['\t'.join(map(str, row)) for row in table_rows]))
        tsv_file.write( '\n')

def save_json(result, output_path):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as w_file:
        json.dump(result, w_file, indent=2)

def save_and_get_input_sequence_text_file_path(input_sequence_fasta_uri, split_inputs_dir=None):
    with urlopen(input_sequence_fasta_uri) as fsrc, tempfile.NamedTemporaryFile(dir=split_inputs_dir, delete=False) as fdst:
        copyfileobj(fsrc, fdst)
        return fdst.name


def commandline_input_prediction(args, parser):
    """ This version takes a file containing an peptide sequences as input."""
    # 1. read input params
    output_prefix = getattr(args, 'output_prefix')
    output_format = getattr(args, 'output_format')

    if output_prefix:
        output_tsv = output_prefix+'.tsv'
        output_json = output_prefix+'.json'
    additional_result_info = {}
    errors = []
    warnings = []
    additional_result_info['warnings'] = warnings

    json_file_name = getattr(args, 'json_filename')
    if json_file_name:
        with open(json_file_name, 'r') as r_file:
            input_data = json.load(r_file)
            validate_result = cluster_validate(input_data)
            errors.extend(validate_result['errors'])
            warnings.extend(validate_result['warnings'])
            peptide_length_range = input_data.get('peptide_length_range', None)
            if peptide_length_range:
                minimum_length, maximum_length = map(int, peptide_length_range)
                lengths = ','.join(map(str,range(minimum_length,maximum_length+1)))
            else:
                minimum_length, maximum_length = 8, 15
                lengths = ''
        if 'input_sequence_text_file_path' in input_data:
            fname = input_data['input_sequence_text_file_path']
            seq_file_type = 'fasta'
            peptide_length_range = input_data['peptide_length_range']
        elif 'input_sequence_fasta_uri' in input_data:
            fname = save_and_get_input_sequence_text_file_path(input_data['input_sequence_fasta_uri'])
            seq_file_type = 'fasta'
            peptide_length_range = input_data['peptide_length_range']
        elif 'input_sequence_text' in input_data:
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_fasta_file:
                fname = tmp_fasta_file.name
                seq_file_type = 'fasta'
                tmp_fasta_file.write(input_data.get('input_sequence_text'))
        elif 'peptide_file_path' in input_data:
            fname = input_data['peptide_file_path']
            seq_file_type = 'peptides'
        else:
            peptide_list = input_data.get('peptide_list')
            seq_file_type = 'peptides'

            #print(peptide_list)
            if not getattr(args, 'assume_valid_flag') and maximum_length:
                for peptide in peptide_list:
                    #print(peptide, len(peptide))
                    if len(peptide) > maximum_length or len(peptide) < minimum_length:
                        peptide_list.remove(peptide)
                        warnings.append('peptide "%s" length is out of valid range (%s)' % (peptide, '%d-%d' % (minimum_length,maximum_length)))
            to_delete = []
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_peptides_file:
                fname = tmp_peptides_file.name
                to_delete.append(fname)
                tmp_peptides_file.write('\n'.join(peptide_list))

        method = input_data.get('method')
        cluster_pct_identity = input_data.get('cluster_pct_identity', 0.97)
        input_allele = input_data.get('alleles')
        additional_result_info['warnings'] = warnings
        additional_result_info["results"] = []

        logging.info('run cluster_analysis: %s, pct_identity: %s, peptide_length_range: %s, method: %s', fname, cluster_pct_identity, peptide_length_range, method)
        result = cluster_analysis(fname, cluster_pct_identity, peptide_length_range, method)

        if output_prefix:
            if output_format.lower()=='tsv':
                result = [result[0]['table_columns']] + result[0]['table_data']
                truncate_file(output_tsv)
                save_tsv(result, output_tsv)
            elif output_format.lower()=='json':
                result[0]['table_columns'] = ["cluster."+ column_name for column_name in result[0]['table_columns']]
                additional_result_info["results"].extend(result)
                save_json(additional_result_info, output_json)
            else:
                eprint('invalid output format: %s' % output_format)
                return

        else:
            if output_format.lower()=='tsv':
                result = [result[0]['table_columns']] + result[0]['table_data']
                print_result(result)
            elif output_format.lower()=='json':
                print(json.dumps(result, indent=2))
            else:
                eprint('invalid output format: %s' % output_format)
                return


def main():
    arg_parser = ClusterArgumentParser()
    args, parser = arg_parser.parse_args()
    arg_parser.aggregate_parameters(args)
    arg_parser.split_parameters(args)

    # cluster prediction
    commandline_input_prediction(args, parser)


if __name__=='__main__':
    main()
