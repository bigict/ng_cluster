# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 20:51:16 2017

@author: bpeters
@modified: sdhanda
"""

import os
import functools
import multiprocessing as mp
import random

num_process = None
if 'epitope_cluster_num_process' in os.environ:
    num_process = int(os.environ['epitope_cluster_num_process'])
chunksize = 1
if 'epitope_cluster_chunksize' in os.environ:
    chunksize = int(os.environ['epitope_cluster_chunksize'])

def _generate_matrix(self, seqs, i, identity_threshold=0.7):
    list_result = [(seqs[i][0], seqs[i][0], float(1))]
    for j in range(i + 1, len(seqs)):
        if seqs[i][1] < seqs[j][1]:
            short_seq, short_len, short_comp = seqs[i]
            long_seq, long_len, long_comp = seqs[j]
        else:
            short_seq, short_len, short_comp = seqs[j]
            long_seq, long_len, long_comp = seqs[i]
        threshold_matches = int((short_len - 0.1) * identity_threshold) + 1
        max_matches = 0
        if self.composition_match(short_comp, long_comp) >= threshold_matches:
            for offset in range(threshold_matches - short_len,
                                long_len - threshold_matches + 1):  # your code is missing similarities that require offsets...
                start_short_pos = max(0, -offset)
                end_short_pos = short_len - max(0, offset + short_len - long_len)
                matches = 0
                for pos in range(start_short_pos, end_short_pos):
                    if short_seq[pos] == long_seq[pos + offset]:
                        matches += 1
                    if not matches + end_short_pos - pos > threshold_matches:
                        break
                else:
                    if matches >= max_matches:
                        max_matches = matches

                list_result.append((seqs[i][0], seqs[j][0], max_matches / float(
                    short_len)))  # I used the same output gathering as you did; this actually takes substantial time. I would only store the connected values, and reserve sufficient space for the list, and/or write to a file
    return list_result

class seq2matrix():
    def composition_match(self, comp1, comp2):
        max_possible_matches = 0
        for c1, c2 in zip(comp1, comp2):
            max_possible_matches += min(c1, c2)
        return max_possible_matches

    def in2seq(self, file):
        seqs = []
        for line in file:
            line = line.strip()
            if line == "":
                continue
            if line[0] == ">":
                continue
            aa_composition = [0] * 20
            for letter in line:
                index = "ACDEFGHIKLMNPQRSTVWY".find(letter)
                if index == -1:
                    print("Invalid AA in sequence", line)
                aa_composition[index] += 1

            seqs.append((line, len(line), aa_composition))
        return seqs

    def generate_matrix(self, seqs, identity_threshold=0.7):

        list_result = []
        with mp.Pool(processes=num_process) as p:
            func = functools.partial(_generate_matrix, self, seqs, identity_threshold=identity_threshold)

            idx_list = list(range(len(seqs)))
            random.shuffle(idx_list)

            for results in p.imap(func, idx_list, chunksize=chunksize):
              list_result += results
        return list_result

        for i in range(len(seqs)):
            list_result.append((seqs[i][0], seqs[i][0], float(1)))
            for j in range(i + 1, len(seqs)):
                if seqs[i][1] < seqs[j][1]:
                    short_seq, short_len, short_comp = seqs[i]
                    long_seq, long_len, long_comp = seqs[j]
                else:
                    short_seq, short_len, short_comp = seqs[j]
                    long_seq, long_len, long_comp = seqs[i]
                threshold_matches = int((short_len - 0.1) * identity_threshold) + 1
                max_matches = 0
                if self.composition_match(short_comp, long_comp) >= threshold_matches:
                    for offset in range(threshold_matches - short_len,
                                        long_len - threshold_matches + 1):  # your code is missing similarities that require offsets...
                        start_short_pos = max(0, -offset)
                        end_short_pos = short_len - max(0, offset + short_len - long_len)
                        matches = 0
                        for pos in range(start_short_pos, end_short_pos):
                            if short_seq[pos] == long_seq[pos + offset]:
                                matches += 1
                            if not matches + end_short_pos - pos > threshold_matches:
                                break
                        else:
                            if matches >= max_matches:
                                max_matches = matches

                        list_result.append((seqs[i][0], seqs[j][0], max_matches / float(
                            short_len)))  # I used the same output gathering as you did; this actually takes substantial time. I would only store the connected values, and reserve sufficient space for the list, and/or write to a file
        return list_result
