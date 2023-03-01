#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--original_alignment', type=str, required=True,
                    help='original_alignment')
parser.add_argument('-o', '--cleaned_alignment', type=str, required=True,
                    help='Out alignment')
args = parser.parse_args()

# ### TESTING VARIABLES
# class Namespace:
#     def __init__(self, **kwargs):
#         self.__dict__.update(kwargs)
# 
# args = Namespace(original_alignment = 'data/extracted_fastas/2023/Novosib_75_realign_trimmed.fasta', cleaned_alignment = 'data/extracted_fastas/2023/cleaned_Novosib_75_realign_trimmed.fasta')
# ### TESTING VARIABLES

alignment_in = AlignIO.read(args.original_alignment, 'fasta')
og_len = len(alignment_in)

# count number of Ds at each position
D_sum = [0] * len(alignment_in[0].seq)
for j in range(og_len):
  for i in range(len(D_sum)):
    if str(alignment_in[j].seq)[i] == "D":
      D_sum[i] += 1
# count number of Es at each position
E_sum = [0] * len(alignment_in[0].seq)
for j in range(og_len):
  for i in range(len(E_sum)):
    if str(alignment_in[j].seq)[i] == "E":
      E_sum[i] += 1
# count number of gaps at each position
gap_sum = [0] * len(alignment_in[0].seq)
for j in range(og_len):
  for i in range(len(gap_sum)):
    if str(alignment_in[j].seq)[i] == "-":
      gap_sum[i] += 1


# find position of 3 most D's/E's
D_max=sorted(D_sum, reverse=True)[:3]
if D_max[2]<max(E_sum):
  D_max[2]=max(E_sum)
  D_switch=False
else:
  D_switch=True

# determine which og positions have over 90% Ds or Es
DDE_pos=[]
for i in range(len(D_sum)):
  if D_sum[i] in D_max:
    DDE_pos.append(i)
  if D_switch is False:
    if E_sum[i] in D_max:
      DDE_pos.append(i)
print(DDE_pos)

# individually check each sequence in alignment to determine if sequences have Ds or Es at correct position
with open(args.cleaned_alignment, 'w') as handle:
  for i in range(len(alignment_in)):
    for j in range(len(DDE_pos)):
      # create list to store if DDE in correct positions
      if j == 0:
        DDE_check = []
      if D_switch is False:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos[j]] == 'D' or str(alignment_in[i].seq)[DDE_pos[j]] == 'E')
      else:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos[j]] == 'D')
    # write to file if all necessary positions have Ds or Es
    if all(DDE_check):
      SeqIO.write(alignment_in[i], handle=handle, format="fasta-2line") # append to file if all is good
