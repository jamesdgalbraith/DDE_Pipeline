#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from re import sub
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--original_alignment', type=str, required=True,
                    help='original_alignment')
parser.add_argument('-o', '--cleaned_alignment', type=str, required=True,
                    help='Out alignment')
args = parser.parse_args()

# ## TESTING VARIABLES
# class Namespace:
#     def __init__(self, **kwargs):
#         self.__dict__.update(kwargs)
# 
# args = Namespace(original_alignment = 'data/extracted_fastas/cdhit_40/renamed_alignments/MuDR.fasta', cleaned_alignment = 'data/extracted_fastas/cdhit_40/DDE_check/MuDR.fasta')
# ## TESTING VARIABLES

## read in data and calculate counters
print(sub("_.*", "", sub(".*/", "", args.original_alignment)))
alignment_in = AlignIO.read(args.original_alignment, 'fasta')
seq_count = len(alignment_in)
aln_length=len(alignment_in[0].seq)

# count number of Ds at each position
D_sum = [0] * aln_length
for j in range(seq_count):
  for i in range(aln_length):
    if str(alignment_in[j].seq)[i] == "D":
      D_sum[i] += 1

# count number of Es at each position
E_sum = [0] * aln_length
for j in range(seq_count):
  for i in range(aln_length):
    if str(alignment_in[j].seq)[i] == "E":
      E_sum[i] += 1

# count number of gaps at each position
nogap_sum = [0] * aln_length
for j in range(seq_count):
  for i in range(aln_length):
    if str(alignment_in[j].seq)[i] != "-":
      nogap_sum[i] += 1

# count number of Ds or Es at each position
DE_sum = [0] * aln_length
for i in range(aln_length):
  DE_sum[i]=D_sum[i]+E_sum[i]

# determine prop of D taking gaps into account
D_sum_nogap=[0] * aln_length
for i in range(aln_length):
  D_sum_nogap[i]=D_sum[i]/nogap_sum[i]

# determine prop of E taking gaps into account
E_sum_nogap=[0] * aln_length
for i in range(aln_length):
  E_sum_nogap[i]=E_sum[i]/nogap_sum[i]

# determine prop of D/E taking gaps into account
DE_sum_nogap=[0] * aln_length
for i in range(aln_length):
  DE_sum_nogap[i]=DE_sum[i]/nogap_sum[i]

# find position of 2 most Ds, ensure over 80% D
D_sum_nogap_max=sorted(D_sum_nogap, reverse=True)[:3]
D_sum_nogap_max=list(filter(lambda score: score >= 0.7, D_sum_nogap_max))
DDE_pos=[]
for i in range(aln_length):
  if D_sum_nogap[i] in D_sum_nogap_max:
    DDE_pos.append(i)
DDE_pos=DDE_pos[0:2]

# set all DEs before second D to be 0
for i in range(DDE_pos[1]+1):
  DE_sum_nogap[i]=0

# find following most common D/E
for i in range(aln_length):
  if DE_sum_nogap[i] == max(DE_sum_nogap):
    DDE_pos.append(i)

print("DDE positions = "+str(DDE_pos))

# individually check each sequence in alignment to determine if sequences have Ds or Es at correct position
with open(args.cleaned_alignment, 'w') as handle:
  for i in range(seq_count):
    for j in range(3):
      # create list to store if DDE in correct positions
      if j == 0:
        DDE_check = []
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos[j]] == 'D')
      if j == 1:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos[j]] == 'D')
      if j == 2:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos[j]] == 'D' or str(alignment_in[i].seq)[DDE_pos[j]] == 'E')
    # write to file if all necessary positions have Ds or Es
    if all(DDE_check):
      SeqIO.write(alignment_in[i], handle=handle, format="fasta-2line") # append to file if all is good

# # Can use for "2 out of 3 ain't bad" approach
# sum(DDE_check)>1
