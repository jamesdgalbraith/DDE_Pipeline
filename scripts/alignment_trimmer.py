#!/usr/bin/env python

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed.')
parser.add_argument('-o', '--out_seq', type=str, required=True,
                    help='Output path for trimmed MSA')
parser.add_argument('-t', '--threshold', type=float, default = 0.5,
                    help='Threshold coverage for inclusion')
args = parser.parse_args()

in_seq = AlignIO.read(args.in_seq, "fasta")

# count bases in each column
for i in range(len(in_seq)):
    for j in range(len(in_seq[0])):
        if j == 0:
            new_list = []
        if str(in_seq[i].seq)[j] == "-":
            new_list.append(0)
        else:
            new_list.append(1)
    if i == 0:
        total_mat = [new_list]
    else:
        total_mat.append(new_list)
total_mat_sum = list(np.sum(total_mat, axis=0))

# determine if columns reach threshold
for i in range(len(total_mat_sum)):
    if i == 0:
        meets_threshold = []
    if total_mat_sum[i]/len(in_seq) >= args.threshold:
        meets_threshold.append(i)

# select columns which reach threshold
for i in range(len(in_seq)):
    if i == 0:
        trimmed_alignment = MultipleSeqAlignment([])
    for j in range(len(meets_threshold)):
        if j == 0:
            new_seq = ""
        k = meets_threshold[j]
        new_seq+=(str(in_seq[i].seq)[k])
    trimmed_alignment.append(SeqRecord(Seq(new_seq), name = in_seq[i].name, id = in_seq[i].id, description = ""))

# save alignment
with open(args.out_seq, "w") as handle:
    AlignIO.write(trimmed_alignment, handle, "fasta")

for i in range(len(trimmed_alignment)):
    if i == 0:
        trimmed_alignment_2 = MultipleSeqAlignment([])
    start_gaps = (str(trimmed_alignment[i].seq[0:10]).count('-'))
    end_gaps = (str(trimmed_alignment[i].seq[(len(trimmed_alignment[0])-10):len(trimmed_alignment[0])]).count('-'))
    prop_gaps = (str(trimmed_alignment[i].seq).count('-'))/len(str(trimmed_alignment[i].seq))
    if start_gaps < 10 and end_gaps < 10 and prop_gaps <= 0.5:
        trimmed_alignment_2.append(SeqRecord(trimmed_alignment[i].seq, name = trimmed_alignment[i].name, id = trimmed_alignment[i].id, description = ""))

index = args.out_seq.rfind('/')
trimmed_out_seq = args.out_seq[:index] + "/extra_" + args.out_seq[index+1:]

# save alignment
with open(trimmed_out_seq, "w") as handle:
    AlignIO.write(trimmed_alignment_2, handle, "fasta")
