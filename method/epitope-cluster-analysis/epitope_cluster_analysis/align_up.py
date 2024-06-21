#!/usr/bin/python
import os

from logging import getLogger

logger = getLogger()




class msa():
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

    def list2fafile(self, peplist, cluster2_tmpdir):
        import tempfile
        if not os.path.exists(cluster2_tmpdir):
            os.makedirs(cluster2_tmpdir)
        filename = tempfile.NamedTemporaryFile(delete=False, suffix=".fa", dir=cluster2_tmpdir).name
        with open(filename, "wt") as tf1:
            seqnum = 0
            for seq in peplist:
                seqnum += 1
                tf1.write(">seq" + str(seqnum) + "\n" + seq + "\n")
        return filename

    def peplist2mat(self, peplist, threshold, cluster2_tmpdir=''):
        from . import identity_matrix_up
        import pandas as pd
        import tempfile

        # set default cluster_tmpdir under MEDIA_ROOT
        if not cluster2_tmpdir:
            tmpdir = os.path.join(MEDIA_ROOT, 'tmp')
            cluster2_tmpdir = os.path.join(tmpdir, 'cluster2')
        # create dir if not exist (on worker servers)
        if not os.path.exists(cluster2_tmpdir):
            os.makedirs(cluster2_tmpdir)

        tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix=".fa", dir=cluster2_tmpdir, mode="w")
        sq_num = 0
        for pep in peplist:
            sq_num += 1
            tmpfile.write(">seq" + str(sq_num) + "\n" + pep + "\n")
        tmpfile.close()
        filename = open(tmpfile.name, 'r')
        c = identity_matrix_up.seq2matrix()
        seqs = c.in2seq(filename)
        out, outcsv = c.generate_matrix(seqs, threshold)
        hdr = ['seq_a', 'seq_b', 'offset', 'shorter', 'longer', 'seqa', 'seqb', 'identity']
        df = pd.DataFrame(outcsv, columns=hdr)
        # df.to_csv(csvfilename,index=False)
        return df

    def update_align(self, peplist, threshold, cluster2_tmpdir):
        import pandas as pd
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt
        import sys
        from operator import itemgetter
        from .identity_matrix_up import seq2matrix
        df = self.peplist2mat(peplist, threshold, cluster2_tmpdir)
        s2m = seq2matrix()
        # df=pd.read_csv(csvfile)
        df = df[df['identity'] >= threshold]
        nxin = df[['seq_a', 'seq_b', 'identity']].values.tolist()
        G = nx.Graph()
        G.add_weighted_edges_from(nxin)
        G1 = nx.Graph([(u, v, d) for (u, v, d) in G.edges(data=True) if d['weight'] >= threshold])
        # sub_graphs = sorted(nx.connected_component_subgraphs(G1), key=len,reverse=True)
        # all_degree = sorted(G1.degree_iter(), key=itemgetter(1), reverse=True)
        all_degree = sorted(dict(G1.degree()).items(), key=itemgetter(1), reverse=True)
        nx.draw(G1, with_labels=True)
        maximal_node = all_degree[0][0]
        cons = maximal_node
        path_information = nx.single_source_shortest_path_length(G1, maximal_node)
        all_paths = set(path_information.values())
        all_paths.remove(0)
        cluster_count = 0
        neighbors = []
        seq1 = []
        seq1.append(maximal_node)
        allseq = []
        for each_path in all_paths:
            next_neighbors = list(k for k, v in path_information.items() if v == each_path)
            cluster_count += 1
            # print each_path,next_neighbors
            next_consensus = []
            ###step 4E.1 append these nodes to make a common consensus
            for each_next_neighbor in next_neighbors:
                # print each_next_neighbor
                # print "cons",cons
                identity_list = s2m.get_identity_pair(cons, each_next_neighbor)
                identity = identity_list[5]
                # print identity_list
                cons1 = identity_list[3]
                pre, post = cons1.split(cons)
                # print "pre",pre,post,cons1,cons
                cons = identity_list[3]
                seq1 = [pre + x + post for x in seq1]
                seq1.append(identity_list[4])
                next_consensus.append(each_next_neighbor)
                neighbors.append(each_next_neighbor)
            fafile = self.list2fafile(seq1, cluster2_tmpdir)
            cons, allseq = self.consensus_seq(fafile)
        pep_order = []
        unorder = [x.replace("-", "") for x in allseq]
        for pep in peplist:
            idx = unorder.index(pep)
            pep_order.append(allseq[idx])
        return cons, pep_order
