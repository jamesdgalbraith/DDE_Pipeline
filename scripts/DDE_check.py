#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--original_alignment', type=str, required=True,
                    help='Existing alignment')
parser.add_argument('-i', '--new_alignment', type=str, required=True,
                    help='New_alignment')
parser.add_argument('-o', '--cleaned_alignment', type=str, required=True,
                    help='Out alignment')
args = parser.parse_args()

print(args.original_alignment)
og_len = len(AlignIO.read(args.original_alignment, 'fasta'))
alignment_in = AlignIO.read(args.new_alignment, 'fasta')
to_check=alignment_in[0:og_len]
to_clean=alignment_in[og_len:len(alignment_in)]

# count number of Ds and Es at each position
DDE_sum = [0] * len(to_check[0].seq)
for j in range(og_len):
  for i in range(len(DDE_sum)):
    if str(to_check[j].seq)[i] == "D" or str(to_check[j].seq)[i] == "E":
      DDE_sum[i] += 1

# determine which og positions have over 90% Ds or Es
DDE_pos=[]
for i in range(len(DDE_sum)):
  if DDE_sum[i]/og_len == 1:
    DDE_pos.append(i)
print(DDE_pos)

# individually check each sequence in alignment to determine if sequences have Ds or Es at correct position
with open(args.cleaned_alignment, 'w') as handle:
  AlignIO.write(to_check, handle, "fasta")
  for i in range(len(to_clean)):
    for j in range(len(DDE_pos)):
      # create list to store if DDE in correct positions
      if j == 0:
        DDE_check = []
      DDE_check.append(str(to_clean[i].seq)[DDE_pos[j]] == 'D' or str(to_clean[i].seq)[DDE_pos[j]] == 'E')
    # write to file if all necessary positions have Ds or Es
    if all(DDE_check):
      SeqIO.write(to_clean[i], handle=handle, format="fasta-2line") # append to file if all is good
