
def format_cluster_table_data(table_data, method):
    # make column 1 a conbined Cluster.Sub-Cluster Number
    # add a new column 'cluster consensus' in the end 
    if method == 'cliques':
        cluster_consensus_dict = {str(row[0]):row[2] for row in table_data if row[1]=='Consensus' or row[1]=='Singleton'}
    else:
        cluster_consensus_dict = {str(row[0])+'.'+str(row[5]):row[2] for row in table_data if row[1]=='Consensus' or row[1]=='Singleton'}

    for row in table_data:
        # method 'cliques' don't have Sub-Cluster Number
        if method == 'cliques':
            if type(row[5]) is list:
                row[5] = ', '.join(row[5])
        else:
            row[0] = str(row[0])+'.'+str(row[5])
            if type(row[6]) is list:
                row[6] = ', '.join(row[6])
            row.pop(5)
        # switch peptide and description columns position to make it same with the header order
        row[-1], row[-2] = row[-2], row[-1]
        row.append(cluster_consensus_dict[str(row[0])])
    #table_data = [list(map(str,row)) for row in table_data]
    return table_data


def cluster_analysis(fasta_file_path, cluster_pct_identity=0.97, peptide_length_range=[0,0], method='cluster-break', peptide_list=None):
    tools_group='cluster'

    threshold = float(cluster_pct_identity)

    if peptide_length_range and len(peptide_length_range) == 2:
        minimum_length, maximum_length = map(int, peptide_length_range)
    else:
        maximum_length = 10000
        minimum_length = 1      
    if not minimum_length:
        minimum_length = 1
    minimum_length = int(minimum_length)      
    if not maximum_length:
        maximum_length = 10000
    maximum_length = int(maximum_length)


    from epitope_cluster_analysis.analysis import process

    result = process(fasta_file_path, threshold, minimum_length, maximum_length, method)
    table_data = result['table_data']
    table_data = format_cluster_table_data(table_data, method)
    visualization_data = result['vis']

    columns = ["cluster_number", "peptide_number", "alignment", "position", "sequence_number", "peptide", "cluster_consensus"]
    if method == 'cliques':
        columns[0] = "clique_number"
    # get final results
    final_results = [
        {
            "type": "peptide_table",
            "table_columns": columns,
            "table_data": table_data
        },
        {
            "type": "network_graph",
            "graph_data": visualization_data
        }
    ]
    
    return final_results
