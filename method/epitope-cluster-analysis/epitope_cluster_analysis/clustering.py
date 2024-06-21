"""
Created on 8 March, 2017

@author: Sandeep Dhanda
"""


class Cluster2():
    import networkx as nx
    import pandas as pd
    import re
    # import matplotlib.pyplot as plt
    # from networkx.drawing.nx_agraph import graphviz_layout
    # import sys
    def consensus_seq(self, filename):
        from Bio import AlignIO
        from Bio.Align import AlignInfo
        alignment = AlignIO.read(open(filename), "fasta")
        all_seq = []
        for record in alignment:
            all_seq.append(str(record.seq))
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = str(summary_align.dumb_consensus(threshold=0.51))
        return consensus, all_seq

    def align(self, peptide_list):
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        from Bio.Seq import Seq
        # from Bio.Alphabet import IUPAC
        from Bio.Align.Applications import ClustalOmegaCommandline
        import tempfile

        in_file = tempfile.NamedTemporaryFile(delete=False, suffix='.fa').name
        out_file = tempfile.NamedTemporaryFile(delete=False, suffix='.aln').name
        counter = 1
        pep_seq = []

        for cn in peptide_list:
            # pep_seq.append(SeqRecord(Seq(cn, IUPAC.protein), id="seq" + str(counter)))
            pep_seq.append(SeqRecord(Seq(cn), id="seq" + str(counter), annotations={'molecule_type': 'protein'}))

            counter += 1
        with open(in_file, "w") as output_handle:
            SeqIO.write(pep_seq, output_handle, "fasta")
        # clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True, force=True)
        clustalomega_cline = ClustalOmegaCommandline("/opt/clustalo/clustalo-1.2.4-Ubuntu-x86_64", infile=in_file,
                                                     outfile=out_file, verbose=True, auto=True, force=True)
        # clustalomega_cline = ClustalOmegaCommandline("/Users/sdhanda/software/clustalw2/clustalw-2.1/src/clustalw2", infile=in_file, outfile=out_file, verbose=True, auto=True, force=True)
        # print clustalomega_cline
        clustalomega_cline()
        cons2, all_seq2 = self.consensus_seq(out_file)
        return cons2, all_seq2

    def get_identity(self, seq_a, seq_b):
        if len(seq_a) > len(seq_b):
            longer = str(len(seq_b) * '-') + seq_a + str(len(seq_b) * '-')
            short_pep = seq_b
        else:
            longer = str(len(seq_a) * '-') + seq_b + str(len(seq_a) * '-')
            short_pep = seq_a
        # print("longer=",longer," short=",short_pep)
        id_list = []
        short_list = []
        for i in range(0, len(longer)):
            longer_pep = longer[i:i + len(short_pep)]
            score = 0;
            for aa, bb in zip(longer_pep, short_pep):
                if aa == bb:
                    score += 1
                # print aa,bb,score
            identity = round(score / float(len(short_pep)), 2)
            short_list.append(longer_pep)
            id_list.append(identity)
        # print("long",longer_pep,"shorter",short,str(identity))
        mx = max(id_list)
        # idx=id_list.index(mx)
        return mx

    def network2consensus(self, G, threshold, cluster2_tmpdir):
        import sys
        import networkx as nx
        from operator import itemgetter
        from .align_up import msa

        # py_ver_major = sys.version_info.major
        py_ver_minor = sys.version_info.minor

        updated = msa()
        # extend=extend_consensus()
        all_consensus = []
        peptide_in_consensus = []
        visualization = []
        # nx.draw(G)
        # plt.show()

        # Assuming that cluster code will only run under 3.8 and above
        if py_ver_minor < 9:
            sub_graphs = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
        else :
            # These codes are for python 3.9+
            G = nx.Graph(G)
            # Find the connected components
            connected_components = list(nx.connected_components(G))

            # Create subgraphs for the connected components
            subgraphs = [G.subgraph(component) for component in connected_components]
            sub_graphs = sorted(subgraphs, key=len, reverse=True)

        for i, sg in enumerate(sub_graphs):
            if len(list(sg.nodes())) > 1:

                # pos = graphviz_layout(sg)
                # all_degree = sorted(sg.degree_iter(), key=itemgetter(1),
                all_degree = sorted(dict(sg.degree()).items(), key=itemgetter(1),
                                    reverse=True)  ## sort the nodes, based on their connectivity
                maximal_node = all_degree[0][0]
                # print all_degree

                #### step 4A find all the neighbors of maximum connectivity node
                neighbors = list(nx.all_neighbors(sg, maximal_node))
                # neighbors.append(maximal_node)
                # print neighbors

                ### step 4B find the alignment and consensus of maximally connected node
                ####cons,neighbor_align=self.align(neighbors)
                cons, neighbor_align = updated.update_align(neighbors, threshold, cluster2_tmpdir)
                # print "consensus",cons
                # print neighbor_align
                ### Step 4 AA createing neighbor list with their edges starts
                neighbor_edge = []
                for each_neighbor in neighbors:
                    connection = (maximal_node, each_neighbor)
                    neighbor_edge.append(connection)
                #### neighbors list ends

                #### step 4C find all the non-neighbors of maximum connectivity node
                # non_neighbors=list(nx.non_neighbors(sg,maximal_node))
                # print "non neighbors", non_neighbors

                ### step 4D get the path information with maximal node
                path_information = nx.single_source_shortest_path_length(sg, maximal_node)

                all_paths = set(path_information.values())
                # print sorted(path_information.items(), key=itemgetter(1))
                # print sorted(path_information), key=path_information.get)
                # print path_information

                # print all_paths
                all_paths.remove(0)  ### remove maximal node sequence
                all_paths.remove(1)  ### remove all first neighbors
                # print cons

                unaligned = []
                ### step 4E follow the nodes accoding to their path
                cluster_count = 1
                for each_path in all_paths:
                    next_neighbors = list(k for k, v in path_information.items() if v == each_path)
                    cluster_count += 1
                    # print each_path,next_neighbors
                    next_consensus = []
                    ######next_consensus.append(cons)
                    # print neighbors
                    # next_consensus.append(neighbors)
                    ###step 4E.1 append these nodes to make a common consensus
                    for each_next_neighbor in next_neighbors:
                        # print each_next_neighbor
                        identity = self.get_identity(cons, each_next_neighbor)
                        if identity >= threshold:
                            next_consensus.append(each_next_neighbor)
                            neighbors.append(each_next_neighbor)
                            # print "peptide added",cons,each_next_neighbor
                        else:
                            unaligned.append(each_next_neighbor)
                            # print "peptide not added", each_next_neighbor
                    # print unaligned
                    if len(next_consensus) > 1:
                        # print next_consensus
                        # print neighbors
                        #### cons,neighbor_align=self.align(neighbors)
                        cons, neighbor_align = updated.update_align(neighbors, threshold, cluster2_tmpdir)
                        # print "next consensus", cons

                # print "aligned",neighbors
                # print cons
                all_consensus.append(cons)
                peptide_in_consensus.append(neighbors)
                # print "unaligned",set(unaligned)

                G.remove_nodes_from(neighbors)
                if not nx.is_frozen(sg):
                  sg.remove_nodes_from(neighbors)
            else:
                single = list(sg.nodes())
                all_consensus.append(single)
                peptide_in_consensus.append(single)
                G.remove_nodes_from(single)

        consensus_peptide = zip(all_consensus, peptide_in_consensus)

        if 9 <= py_ver_minor:
            # Freeze the graph for python version 3.9+
            G = nx.freeze(G)

        return G, consensus_peptide, visualization

    def graph2output(self, G, threshold, cluster_breaking_identity, cluster2_tmpdir):
        import sys
        import networkx as nx
        import pandas as pd
        import re
        from .align_up import msa
        updated = msa()

        # py_ver_major = sys.version_info.major
        py_ver_minor = sys.version_info.minor

        if py_ver_minor < 9 :
            sub_graphs = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)  ## First sub-network has the highest connectivity
        else :
            G = nx.Graph(G)
            # Find the connected components
            connected_components = list(nx.connected_components(G))

            # Create subgraphs for the connected components
            subgraphs = [G.subgraph(component) for component in connected_components]
            sub_graphs = sorted(subgraphs, key=len, reverse=True)

        nested_list = []
        vis = []
        clustering_counter = []
        counter1 = 0
        for i, subgraph in enumerate(sub_graphs):
            counter1 += 1
            counter2 = 0

            while len(subgraph.nodes()) > 0:
                # print len(subgraph.nodes())
                # nx.draw(G,with_labels=True)
                # plt.show()

                subgraph, consensus_peptides, visualization = self.network2consensus(subgraph,
                                                                                     cluster_breaking_identity,
                                                                                     cluster2_tmpdir)
                # vis.append(visualization)
                for cons, peptides in consensus_peptides:
                    counter2 += 1
                    counter = float(str(counter1) + "." + str(counter2))
                    cluster_count2 = str(counter1) + "_" + str(counter2)
                    clustering_counter.append(cluster_count2)
                    all_peptide = []
                    if len(peptides) > 1:
                        ##cons2,neighbor_align=self.align(peptides)
                        cons2, neighbor_align = updated.update_align(peptides, threshold, cluster2_tmpdir)
                        cons_name = "Consensus"
                        align_name = "-"
                        start_pos = 0
                    else:
                        cons2 = peptides[0]
                        neighbor_align = peptides
                        cons_name = "Singleton"
                        align_name = peptides[0]
                        start_pos = '-'
                    edge_list = G.edges(nbunch=peptides, data=True)
                    for edge in edge_list:
                        if edge[0] and edge[1] in peptides:
                            rs_dict = {}
                            rs_dict["source"] = edge[0]
                            rs_dict["target"] = edge[1]
                            rs_dict["cluster"] = [cluster_count2]
                            # res_var="{source : \"" + edge[0] + "\", target : \"" + edge[1] + "\", weight : \""+str(edge[2].values()[0])+"\", cluster : \"" +str(cluster_count2)+"\"},\n"
                            # vis.append(res_var)
                            vis.append(rs_dict)

                    consensus = [counter1, counter2, cons_name, str(cons2), start_pos, '-', align_name]
                    nested_list.append(consensus)
                    if len(peptides) > 1:
                        # all_nodes=G.nodes()

                        # res = list(set(all_nodes)^set(peptides))
                        # G2=G.copy()
                        # G2.remove_nodes_from(res)
                        # print G2.nodes()

                        for a, peptide in enumerate(peptides):
                            all_peptide.append(peptide)

                            for ret in re.finditer('\w+', neighbor_align[a]):
                                start = ret.start() + 1
                            # res=str(counter)+","+str(cons2)+","+str(len(peptides))+ "," +str(neighbor_align[a])+","+str(start)+","+str(len(peptide))+","+peptide+"\n"
                            listed_result = [counter1, counter2, a + 1, neighbor_align[a], start, len(peptide), peptide]
                            # print res
                            nested_list.append(listed_result)
                            # print listed_result
                            # ft.write(res)

                        ##### to plot all the nodes in subgraph   Uncomment next four lines
            # sg=H.subgraph(all_peptide)
            # print sg.edges()
            # nx.draw(sg,with_labels=True)
            # plt.show()

        result = pd.DataFrame(nested_list)
        # print result
        result.sort_values([0, 1, 4, 5], ascending=[True, True, True, False], inplace=True)
        result.columns = ['Cluster_number', 'Sub_cluster_number', 'Consensus', 'Alignment', 'Position', 'Length',
                          'Peptide']
        result['Peptide_number'] = result.groupby(['Cluster_number', 'Sub_cluster_number'])['Cluster_number'].transform(
            lambda x: pd.Series(range(0, len(x)), index=x.index))

        # same with 0.24.2 version code: result.Peptide_number.replace(0, result.Consensus, inplace=True)
        result['Peptide_number'] = result['Peptide_number'].astype('object')
        for index, row in result.iterrows():
            if row['Peptide_number'] == 0:
                result.at[index, 'Peptide_number'] = row['Consensus']

        # result['Peptide_number'].loc[result.Consensus == 0] = "-"
        result = result[['Cluster_number', 'Peptide_number', 'Alignment', 'Position', 'Peptide', 'Sub_cluster_number']]
        return result, vis, clustering_counter

    def cliques2output(self, G, threshold, cluster_breaking_identity, cluster2_tmpdir):
        import networkx as nx
        import pandas as pd
        import re
        from .align_up import msa
        updated = msa()
        cliques_listed = sorted(list(nx.find_cliques(G)), key=len,
                                reverse=True)  ## First sub-network has the highest connectivity
        nested_list = []
        vis = []
        done = {}
        res_dict2 = {}
        clustering_count = []
        counter1 = 0
        # with open(final_results,"at") as ft:
        for i, sublist in enumerate(cliques_listed):
            counter1 += 1
            counter2 = 0
            # print sublist
            if len(sublist) > 0:
                # print len(subgraph.nodes())
                # nx.draw(G,with_labels=True)
                # plt.show()
                # subgraph=nx.from_edgelist(sublist)
                sblist = []
                for a in sublist:
                    for b in sublist:
                        ab = (a, b, 1)
                        sblist.append(ab)
                subgraph = nx.Graph()
                subgraph.add_weighted_edges_from(sblist)
                # print subgraph.nodes()
                subgraph, consensus_peptides, visualization = self.network2consensus(subgraph,
                                                                                     cluster_breaking_identity,
                                                                                     cluster2_tmpdir)
                # vis.append(visualization)
                for cons, peptides in consensus_peptides:
                    counter2 += 1
                    counter = float(str(counter1) + "." + str(counter2))
                    cluster_count2 = str(counter1) + "_"
                    clustering_count.append(cluster_count2)
                    all_peptide = []
                    if len(peptides) > 1:
                        ####cons2,neighbor_align=self.align(peptides)
                        cons2, neighbor_align = updated.update_align(peptides, threshold, cluster2_tmpdir)
                        cons_name = "Consensus"
                        align_name = "-"
                        start_pos = 0
                    else:
                        cons2 = peptides[0]
                        neighbor_align = peptides
                        cons_name = "Singleton"
                        align_name = peptides[0]
                        start_pos = '-'
                    edge_list = G.edges(nbunch=peptides, data=True)
                    for edge in edge_list:
                        if edge[0] and edge[1] in peptides:
                            # res_var="{source : \"" + edge[0] +"_" + str(counter1) + "\", target : \"" + edge[1] + "_" +str(counter1)  + "\", weight : \""+str(edge[2].values()[0])+"\", cluster : \"" +str(cluster_count2)+"\"},\n"
                            res_var = "{source : \"" + edge[0] + "\", target : \"" + edge[1] + "\", weight : \"" + str(
                                list(edge[2].values())[0]) + "\", cluster : \"" + str(cluster_count2) + "\"}"
                            res_dict = {}

                            res_dict2.setdefault(edge[0] + "_" + edge[1], []).append(cluster_count2)
                            res_dict2.setdefault(edge[1] + "_" + edge[0], []).append(cluster_count2)
                        # res_var="{source : \"" + edge[0] + "\", target : \"" + edge[1]  + "\", weight : \""+str(edge[2].values()[0])+"\", cluster : \"" +str(cluster_count2)+"\"}"
                        # vis.append(res_var)
                        # vis.append(res_dict)

                    consensus = [counter1, cons_name, str(cons2), start_pos, '-', align_name]
                    nested_list.append(consensus)
                    if len(peptides) > 1:
                        for a, peptide in enumerate(peptides):
                            all_peptide.append(peptide)

                            for ret in re.finditer('\w+', neighbor_align[a]):
                                start = ret.start() + 1
                            # res=str(counter)+","+str(cons2)+","+str(len(peptides))+ "," +str(neighbor_align[a])+","+str(start)+","+str(len(peptide))+","+peptide+"\n"
                            listed_result = [counter1, a + 1, neighbor_align[a], start, len(peptide), peptide]
                            # print res
                            nested_list.append(listed_result)
                            # print listed_result
                            # ft.write(res)

                        ##### to plot all the nodes in subgraph   Uncomment next four lines
                        # sg=H.subgraph(all_peptide)
                        # print sg.edges()
                        # nx.draw(sg,with_labels=True)
                        # plt.show()
        for k, v in res_dict2.items():
            rs_dict = {}
            rs_dict["source"] = k.split("_")[0]
            rs_dict["target"] = k.split("_")[1]
            rs_dict["cluster"] = v
            vis.append(rs_dict)
        result = pd.DataFrame(nested_list)
        # print vis
        result.sort_values([0, 3, 4], ascending=[True, True, False], inplace=True)
        result.columns = ['Cluster_number', 'Consensus', 'Alignment', 'Position', 'Length', 'Peptide']
        result['Peptide_number'] = result.groupby('Cluster_number')['Cluster_number'].transform(
            lambda x: pd.Series(range(0, len(x)), index=x.index))
        # same with 0.24.2 version code: result.Peptide_number.replace(0, result.Consensus, inplace=True)
        result['Peptide_number'] = result['Peptide_number'].astype('object')
        for index, row in result.iterrows():
            if row['Peptide_number'] == 0:
                result.at[index, 'Peptide_number'] = row['Consensus']

        # result['Peptide_number'].loc[result.Consensus == 0] = "-"
        result = result[['Cluster_number', 'Peptide_number', 'Alignment', 'Position', 'Peptide']]
        return result, vis, clustering_count
