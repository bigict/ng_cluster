# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 20:51:16 2017

@author: bpeters
@modified: sdhanda
"""


class seq2matrix():
    def composition_match(self, comp1, comp2):
        max_possible_matches = 0
        for c1, c2 in zip(comp1, comp2):
            max_possible_matches += min(c1, c2)
        return max_possible_matches

    def in2seq(self, file):
        seqs = []
        for line in file:
            line = line.strip().upper()
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
        list_result_csv = []
        for i in range(len(seqs)):
            # list_result.append((seqs[i][0], seqs[i][0],float(1)))
            seq_i = seqs[i][0]
            seq_j = seqs[i][0]
            list_result.append((seqs[i][0], seqs[i][0], float(1)))
            list_result_csv.append((seqs[i][0], seqs[i][0], 0, seqs[i][0], seqs[i][0], seq_i, seq_j, float(1)))
            for j in range(i + 1, len(seqs)):
                seq_j = seqs[j][0]
                if seqs[i][1] < seqs[j][1]:
                    short_seq, short_len, short_comp = seqs[i]
                    long_seq, long_len, long_comp = seqs[j]
                else:
                    short_seq, short_len, short_comp = seqs[j]
                    long_seq, long_len, long_comp = seqs[i]
                threshold_matches = int((short_len - 0.1) * identity_threshold) + 1
                max_matches = 0
                max_offset = 0
                seqs_a = short_seq
                seqs_b = long_seq
                if self.composition_match(short_comp, long_comp) >= threshold_matches:
                    for offset in range(threshold_matches - short_len, long_len - threshold_matches + 1):
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
                                max_offset = offset
                                # print "offset",offset,short_seq,long_seq
                                if max_offset >= 0:
                                    seqs_a = "-" * max_offset + short_seq
                                    seqs_b = long_seq
                                    if long_seq == seqs[i][0]:
                                        seq_i = seqs_b
                                        seq_j = seqs_a
                                    else:
                                        seq_i = seqs_a
                                        seq_j = seqs_b
                                else:
                                    max_offset1 = abs(max_offset)
                                    seqs_b = "-" * max_offset1 + long_seq
                                    seqs_a = short_seq
                                    if long_seq == seqs[i][0]:
                                        seq_i = seqs_b
                                        seq_j = seqs_a
                                    else:
                                        seq_i = seqs_a
                                        seq_j = seqs_b
                        # list_result.append((seqs[i][0], seqs[j][0],max_matches/float(short_len)))
                list_result.append((seqs[i][0], seqs[j][0], max_matches / float(short_len)))
                list_result_csv.append(
                    (seqs[i][0], seqs[j][0], max_offset, seqs_a, seqs_b, seq_i, seq_j, max_matches / float(short_len)))
                list_result_csv.append((seqs[j][0], seqs[i][0], max_offset * -1, seqs_b, seqs_a, seq_j, seq_i,
                                        max_matches / float(short_len)))
        return list_result, list_result_csv

    def get_identity_pair(self, seq_a, seq_b):
        if len(seq_a) > len(seq_b):
            longer = str(len(seq_b) * '-') + seq_a + str(len(seq_b) * '-')
            short_pep = seq_b
        else:
            longer = str(len(seq_a) * '-') + seq_b + str(len(seq_a) * '-')
            short_pep = seq_a
        id_list = []
        short_list = []
        for i in range(0, len(longer)):
            longer_pep = longer[i:i + len(short_pep)]
            score = 0;
            for aa, bb in zip(longer_pep, short_pep):
                if aa == bb:
                    score += 1
            identity = round(score / float(len(short_pep)), 2)
            short_list.append(longer_pep)
            id_list.append(identity)
        mx = max(id_list)
        idx = id_list.index(mx)
        offset = longer.find(short_list[idx]) - len(short_pep)
        if short_pep == seq_a:
            if offset >= 0:
                ## tested correct
                post_short = len(seq_b) - offset - len(short_pep)
                new_short = "-" * offset + short_pep + "-" * post_short
                post = len(new_short) - len(seq_b)
                new_long = seq_b + "-" * post
                res = [short_pep, seq_b, offset, new_short, new_long, mx]
            else:
                offset1 = -1 * offset
                ## tested correct
                ##print "test short a"
                new_long = "-" * offset1 + seq_b
                post = len(new_long) - len(seq_a)
                new_short = short_pep + "-" * post
                res = [short_pep, seq_b, offset, new_short, new_long, mx]
        # print(short_pep,clean_longer,offset,short_pep,longer,new_short,new_long,mx)
        else:
            if offset >= 0:
                ##tested correct
                # print "test short b"
                post_short = len(seq_a) - offset - len(short_pep)
                new_short = "-" * offset + short_pep + "-" * post_short
                post = len(new_short) - len(seq_a)
                new_long = seq_a + "-" * post
                # res=[seq_a,short_pep,offset,short_pep,longer,new_long,new_short,mx]
                res = [seq_a, short_pep, offset, new_long, new_short, mx]
            else:
                ##print "Not tested short b, verify first"
                offset1 = -1 * offset
                new_long = "-" * offset1 + seq_a
                post = len(new_long) - len(seq_b)
                new_short = short_pep + "-" * post
                res = [seq_a, short_pep, offset, new_long, new_short, mx]
        # print(clean_longer,short_pep,offset,short_pep,longer,new_long,new_short,mx)
        return res
