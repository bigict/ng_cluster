import logging

import networkx as nx

from . import identity_matrix, clustering

def processing(seqs, num, uniq_num, threshold, method, minimum_len, maximum_len, cluster_breaking_identity,
               cluster2_tmpdir, fasta_seq):

    c = identity_matrix.seq2matrix()

    logging.info('epitope_cluster_analysis.processing generate_matrix ....')
    out = c.generate_matrix(seqs, threshold)
    # print out

    logging.info('epitope_cluster_analysis.processing build_graph, num_edges=%d', len(out))
    G = nx.Graph()  # or DiGraph, MultiGraph, MultiDiGraph, etc
    G.add_weighted_edges_from(out)
    del out

    logging.info('epitope_cluster_analysis.processing filter_graph ...')
    G1 = nx.Graph([(u, v, d) for (u, v, d) in G.edges(data=True) if d['weight'] >= threshold])

    logging.info('epitope_cluster_analysis.processing clustering ....')
    ### step 4. find the consensus for nodes
    clusters = clustering.Cluster2()
    if method == 'cluster':
        cluster_breaking_identity = 0
    if method == 'cliques':
        result, visualization, all_clusters = clusters.cliques2output(G1, threshold, cluster_breaking_identity,
                                                                      cluster2_tmpdir)
    else:
        result, visualization, all_clusters = clusters.graph2output(G1, threshold, cluster_breaking_identity,
                                                                    cluster2_tmpdir)
    result['description'] = result['Peptide'].map(fasta_seq)
    # result=result.drop_duplicates('Peptide','Position')#['description']=result['Peptide'].map(fasta_seq)
    # print result
    result.fillna('-', inplace=True)
    result['Position'].loc[result.Position == 0] = '-'
    result_to_django = result.values.tolist()
    result.set_index('Cluster_number', inplace=True)
    # print len(result[result.'Peptide number' == 'Singleton'])
    # print cluster_count
    # print visualization
    # print type(visualization)
    # print all_clusters
    # return render(request,'cluster2/results.html',{'num':num,'uniq':uniq_num,'threshold':threshold*100,'table_data':result_to_django,'vis':visualization,'all_clusters' : all_clusters,'method':method})#,'csvname':csvfile.name})

    result_data = {'num': num, 'uniq': uniq_num, 'threshold': threshold * 100, 'table_data': result_to_django,
                   'vis': visualization, 'all_clusters': all_clusters, 'method': method, 'min_len': minimum_len,
                   'max_len': maximum_len}

    return result_data
