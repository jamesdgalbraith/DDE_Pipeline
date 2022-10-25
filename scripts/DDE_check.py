#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import os
from re import sub
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--original_alignment', type=str, required=True,
                    help='Existing alignment')
parser.add_argument('-i', '--new_alignment', type=str, required=True,
                    help='New_alignment')
parser.add_argument('-o', '--cleaned_alignment', type=str, required=True,
                    help='Out alignment')
args = parser.parse_args()

# original_alignment="starting_seq/Transib.fasta"
# new_alignment="seq_pass_3/Transib.fasta"
# cleaned_alignment="seq_pass_4/Transib.fasta"

og_len = len(AlignIO.read(args.original_alignment, 'fasta'))
alignment_in = AlignIO.read(args.new_alignment, 'fasta')
consensus_seq = str(AlignInfo.SummaryInfo(alignment_in[0:og_len]).gap_consensus(threshold=0.9, ambiguous="X"))

# find positions of Ds and Es
def DDE_find(Cs):
  for i in range(len(Cs)):
    if i == 0:
      DDE_pos=[]
    if Cs[i] == "D" or Cs[i] == "E":
      DDE_pos.append(i)
  return DDE_pos

# Check if seq is correct at each position, input is positions containing 90% Ds&Es and alignment_in[i].seq
def DDE_check(Dp, Al, handle):
  for j in range(len(Dp)):
    if j == 0:
      Dc=[] # make list for storing if D or E
    Dc.append(str(Al.seq)[Dp[j]] == consensus_seq[Dp[j]])
    if j == len(Dp)-1 and all(Dc) is True:
      SeqIO.write(Al, handle=handle, format="fasta-2line") # append to file if all is good
print(DDE_find(consensus_seq))
with open(args.cleaned_alignment, 'w') as handle:
  for k in range(len(alignment_in)):
    DDE_check(DDE_find(consensus_seq), alignment_in[k], handle)



  
