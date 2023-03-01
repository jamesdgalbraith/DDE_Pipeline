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

# ### TESTING VARIABLES
# class Namespace:
#     def __init__(self, **kwargs):
#         self.__dict__.update(kwargs)
# 
# args = Namespace(original_alignment = 'data/extracted_fastas/2023_04_threshold/Sola_realigned_trimmed.fasta', cleaned_alignment = 'data/extracted_fastas/2023_04_threshold/Sola_realigned_trimmed_DDE.fasta')
# ### TESTING VARIABLES
print(sub("_.*", "", sub(".*/", "", args.original_alignment)))
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
nogap_sum = [0] * len(alignment_in[0].seq)
for j in range(og_len):
  for i in range(len(nogap_sum)):
    if str(alignment_in[j].seq)[i] != "-":
      nogap_sum[i] += 1

# count number of Ds or Es at each position
DE_sum = [0] * len(alignment_in[0].seq)
for j in range(og_len):
  for i in range(len(DE_sum)):
    if str(alignment_in[j].seq)[i] == "D" or str(alignment_in[j].seq)[i] == "E":
      DE_sum[i] += 1

# find position of 3 most D's/E's
D_max=sorted(D_sum, reverse=True)[:3]
if D_max[2]<max(E_sum):
  D_max[2]=max(E_sum)
  D_switch=False
else:
  D_switch=True

DE_max=sorted(DE_sum, reverse=True)[:3]

# determine which og positions have over 90% Ds or Es
DDE_pos=[]
for i in range(len(D_sum)):
  if D_sum[i] in D_max:
    DDE_pos.append(i)
  if D_switch is False:
    if E_sum[i] in D_max:
      DDE_pos.append(i)
print("DDE pos 1 = "+str(DDE_pos))

DDE_pos_2=[]
for i in range(len(DE_sum)):
  if DE_sum[i] in DE_max:
    DDE_pos_2.append(i)
print("DDE pos 2 = "+str(DDE_pos_2))

# determine prop of D/E taking gaps into account
DE_sum_nogap=[0] * len(alignment_in[0].seq)
for i in range(len(DE_sum)):
  DE_sum_nogap[i]=DE_sum[i]/nogap_sum[i]

# determine max D/E taking gap into account
DE_sum_nogap_max=sorted(DE_sum_nogap, reverse=True)[:3]

# determine position of max D/E taking gap into account
DDE_pos_3=[]
for i in range(len(DE_sum_nogap)):
  if DE_sum_nogap[i] in DE_sum_nogap_max:
    DDE_pos_3.append(i)
print("DDE pos 3 = "+str(DDE_pos_3))

# # manual for Academ
# DDE_pos_3=[5, 171, 410]
# # manual for Sola
# DDE_pos_3=[5, 171, 410]

# individually check each sequence in alignment to determine if sequences have Ds or Es at correct position

with open(sub("\.fasta", "1.fasta", args.cleaned_alignment), 'w') as handle:
  for i in range(len(alignment_in)):
    for j in range(len(DDE_pos)):
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

with open(sub("\.fasta", "2.fasta", args.cleaned_alignment), 'w') as handle:
  for i in range(len(alignment_in)):
    for j in range(len(DDE_pos_2)):
      # create list to store if DDE in correct positions
      if j == 0:
        DDE_check = []
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_2[j]] == 'D')
      if j == 1:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_2[j]] == 'D')
      if j == 2:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_2[j]] == 'D' or str(alignment_in[i].seq)[DDE_pos_2[j]] == 'E')
    # write to file if all necessary positions have Ds or Es
    if all(DDE_check):
      SeqIO.write(alignment_in[i], handle=handle, format="fasta-2line") # append to file if all is good

with open(sub("\.fasta", "3.fasta", args.cleaned_alignment), 'w') as handle:
  for i in range(len(alignment_in)):
    for j in range(len(DDE_pos_3)):
      # create list to store if DDE in correct positions
      if j == 0:
        DDE_check = []
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_3[j]] == 'D')
      if j == 1:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_3[j]] == 'D')
      if j == 2:
        DDE_check.append(str(alignment_in[i].seq)[DDE_pos_3[j]] == 'D' or str(alignment_in[i].seq)[DDE_pos_3[j]] == 'E')
    # write to file if all necessary positions have Ds or Es
    if all(DDE_check):
      SeqIO.write(alignment_in[i], handle=handle, format="fasta-2line") # append to file if all is good
