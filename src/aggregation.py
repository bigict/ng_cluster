#!/usr/bin/env python
# coding: utf-8

import os
import json
import re
import pandas as pd

# read in the job description file, which is where we'll pull the output locations from
# note the '-mod' suffix indicated that I needed to modify the job descriptions file so that the expected
# output paths matched my local system

def save_json(result, output_path):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as w_file:
        json.dump(result, w_file, indent=2)

def read_json_file(json_file_path):
    with open(json_file_path, 'r') as r_file:
        return json.load(r_file)

# get unique_vals for string type columns, and field_ranges for number type for a table
def get_filter_ranges(table):
    warnings = []
    object_indices = []
    numeric_indices = []
    for i,dtype in enumerate(table.dtypes):
        if dtype == "object":
            object_indices.append(i)
        else:
            numeric_indices.append(i)

    all_cols = list(table)
    unique_vals = dict()
    for i in object_indices:
        unique_vals[all_cols[i]] = list(map(lambda x: re.sub("^nan$", "-", str(x)), table.iloc[:,i].unique()))

    field_ranges = dict()
    for i in numeric_indices:
        min_value = str(table.iloc[:,i].min())
        max_value = str(table.iloc[:,i].max())
        f = all_cols[i]
        if min_value != 'nan' and max_value != 'nan':
            field_ranges.setdefault(f, {})
            field_ranges[all_cols[i]]['min'] = min_value
            field_ranges[all_cols[i]]['max'] = max_value
        else:
            warnings.append("can not calculate min & max for %s" % f)

    return unique_vals, field_ranges, warnings

def aggregate_result_file(job_desc_file, output_dir, aggregate_output_dir, has_consensus=None):
    job_desc = read_json_file(job_desc_file)

    all_results = dict()
    warnings = set()
    errors = set()

    for j in job_desc:
        # skip over non-prediction jobs
        if j['job_type'] != 'prediction':
            continue

        # loop through the expected outputs of prediction jobs
        for o in j['expected_outputs']:
            output = read_json_file(o)
            for r in output['results']:
                t = r['type']
                # TODO: deal with non-conforming results here
                if t == 'network_graph':
                    if t not in all_results:
                        all_results[t] = r
                    else:
                        if "graph_data" in all_results[t]:
                            all_results[t]["graph_data"].extend(r["graph_data"])
                else:
                    df = pd.DataFrame(r['table_data'], columns=r['table_columns'])
                    if t not in all_results:
                        # if the prediction type has not yet been included, we need to initialize it
                        all_results[t] = df
                    else:
                        # if we get here, there are already results of this type for this method,
                        # so we concatenate the df
                        all_results[t] = pd.concat([all_results[t], df])

            # keep all unique warnings & errors
            if 'warnings' in output:
                for w in output['warnings']:
                    warnings.add(w)

            if 'errors' in output:
                for e in output['errors']:
                    errors.add(e)

    # Finally, we put everything back together.  Note that the keys are slightly different than what we had previously - 'columnns' instead of 'table_columns' and 'data' instead of 'table_data'.  We can update the names, but we would have to create new variables which seems unnecessarily expensive.  It's probably best to simply update the downstream code that needs to work with this data.

    # now let's put this all into an aggregated object that we can dump as JSON
    final_results = list()
    for t in all_results.keys():
        if t != 'network_graph':
            # drop all na and all "-" columns
            if t in ['peptide_table', 'residue_table']:
                df_t = (all_results[t].fillna('-') == '-').all()
                empty_columns = df_t[df_t]
                if not empty_columns.empty:
                    warnings.add("The following fields in %s were empty and have been removed from the output: %s" % (t,', '.join(empty_columns.keys())))
                all_results[t] =  all_results[t].drop(empty_columns.index, axis=1)

            # get unique_vals for string type columns, and field_ranges for number type for each table
            unique_vals, field_ranges, filter_warnings = get_filter_ranges(all_results[t])
            warnings.update(filter_warnings)
            # replace na with "-" after get filter information (e.g. min, max)
            if t in ['peptide_table', 'residue_table']:
                all_results[t] = all_results[t].fillna('-')
            dict_result = all_results[t].to_dict(orient='split')
            del dict_result['index']
            dict_result['result_type'] = t
            dict_result['table_columns'] = dict_result.pop('columns')
            dict_result['table_data'] = dict_result.pop('data')
            dict_result['unique_vals'] = unique_vals
            dict_result['field_ranges'] = field_ranges
            final_results.append(dict_result)

    # add the processing plots onto the results
    if 'network_graph' in all_results:
        if 'type' in all_results['network_graph']:
            all_results['network_graph']['result_type'] = all_results['network_graph'].pop('type')
        final_results.append(all_results['network_graph'])

    aggregated_results = {
        'warnings': list(warnings),
        'errors': list(errors),
        'results': final_results
    }

    # this is the final aggregated result object that can be written to JSON
    # That's it!  The above object can now be dumped to JSON.

    # save result to output path
    # print('aggregated_results: %s' % aggregated_results)
    final_result_file_path = os.path.join(aggregate_output_dir, 'aggregated_result.json')
    save_json(aggregated_results, final_result_file_path)
    print('aggregated_result_path:%s' % final_result_file_path)

    return

