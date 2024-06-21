
import csv, os
import operator
import tempfile
import logging


from .tasks import processing
from . import identity_matrix


def process(seqfile_name, threshold=0.7, minimum_len=1, maximum_len=100, method='cluster-break'):

    cluster2_tmpdir = tempfile.TemporaryDirectory()
    tmpfile2 = open(cluster2_tmpdir.name + "/processed.fa", 'wt')

    threshold = float(threshold)
    # print(minimum_len)
    minimum_len = int(minimum_len)
    maximum_len = int(maximum_len)

    cluster_breaking_identity = threshold
    fasta = open(seqfile_name, 'r')
    from Bio import SeqIO
    peplist = []
    fasta_seq = {}
    for ids, entry in enumerate(SeqIO.parse(fasta, 'fasta')):
        if ((len(str(entry.seq)) >= minimum_len) and (len(str(entry.seq)) <= maximum_len)):
            peplist.insert(ids, str(entry.seq))
            # fasta_seq[str(entry.seq)]=entry.description
            fasta_seq.setdefault(str(entry.seq), []).append(entry.description)
            tmpfile2.write('>' + str(entry.description) + '\n' + str(entry.seq) + '\n');
    fasta.close()
    tmpfile2.close()
    num = len(peplist)
    uniq_num = len(set(peplist))

    c = identity_matrix.seq2matrix()
    filename = open(cluster2_tmpdir.name + "/processed.fa", 'r')
    seqs = c.in2seq(filename)
    filename.close()

    logging.info('epitope_cluster_analysis.analysis.processing ...')
    result_data = processing(seqs, num, uniq_num, threshold, method, minimum_len, maximum_len,
                                 cluster_breaking_identity, cluster2_tmpdir.name, fasta_seq)

    cluster2_tmpdir.cleanup()
    return result_data


