import json
import os
import shutil
import random
import string

def generate_random_str(length):
    return ''.join(random.sample(string.digits+string.ascii_letters, length))

def save_json(result, output_path):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as w_file:
        json.dump(result, w_file, indent=2)

def split_parameters_file(json_filename, parameters_output_dir=None, split_inputs_dir=None, assume_valid=False):
    
    with open(json_filename, 'r') as r_file:
        input_data = json.load(r_file)

    # recreate the directory
    if os.path.exists(parameters_output_dir):
        shutil.rmtree(parameters_output_dir)

    # if not given, create a folder for it
    if not parameters_output_dir:
        parameters_output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'job_%s/splitted_parameters' % generate_random_str(6)))
    # if not given, use parameters_output_dir for generated sequence files as well
    if not split_inputs_dir:
        split_inputs_dir = parameters_output_dir
    # create dir if not exist:
    os.makedirs(parameters_output_dir, exist_ok=True)
    os.makedirs(split_inputs_dir, exist_ok=True)

    # no need to really split input for cluster analysis
    output_data  = [input_data]

    parameters_output_dir = os.path.abspath(parameters_output_dir) 
    result_output_dir = os.path.abspath(os.path.join(parameters_output_dir, os.pardir, 'results'))
    aggregate_dir = os.path.abspath(os.path.join(parameters_output_dir, os.pardir, 'aggregate'))
    base_dir = os.path.abspath(os.path.join(parameters_output_dir, os.pardir))
    job_descriptions_path = os.path.abspath(os.path.join(base_dir, 'job_descriptions.json'))
    cluster_predict_executable_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'run_cluster.py'))

    os.makedirs(result_output_dir, exist_ok=True)
    job_descriptions = []
    aggregate_depends_on_job_ids = []
    job_id = -1
    for i, data_unit in enumerate(output_data):
        job_id = i
        data_unit_file_path = os.path.join(parameters_output_dir, '%d.json' % i)
        save_json(data_unit, data_unit_file_path)
        shell_cmd='%s -j %s/%s.json -o %s/%s -f json' % (cluster_predict_executable_path, parameters_output_dir, job_id, result_output_dir, job_id)
        if assume_valid:
            shell_cmd += ' --assume-valid'
        job_description = dict(
            shell_cmd=shell_cmd,
            job_id=job_id,
            job_type="prediction",
            depends_on_job_ids=[],
            expected_outputs=['%s/%s.json' % (result_output_dir, job_id)]
        )
        aggregate_depends_on_job_ids.append(job_id)
        job_descriptions.append(job_description)
    print('parameters_output_dir: %s' % os.path.abspath(parameters_output_dir))

    # add aggreate job
    # job_id == -1 means no job is required to run
    if job_id > -1:
        job_id +=1
        shell_cmd='%s --aggregate --job-desc-file=%s --aggregate-input-dir=%s --aggregate-result-dir=%s' % (cluster_predict_executable_path, job_descriptions_path, result_output_dir, aggregate_dir)

        aggreate_job_description = dict(
            shell_cmd=shell_cmd,
            job_id=job_id,
            job_type="aggregate",
            depends_on_job_ids=aggregate_depends_on_job_ids,
            expected_outputs=['%s/aggregated_result.json' % aggregate_dir], 
        )
        job_descriptions.append(aggreate_job_description)
    save_json(job_descriptions, job_descriptions_path)
    print('job_descriptions_path: %s' % os.path.abspath(job_descriptions_path))

    return
