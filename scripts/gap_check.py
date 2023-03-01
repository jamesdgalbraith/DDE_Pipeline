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
parser.add_argument('-t', '--threshold', type=int, default=80,
                    help='Maximum % of sequences allowed to be gaps')                    
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

# write sequences to file if criteria met
with open(args.cleaned_alignment, 'w') as handle:
  for i in range(seq_count):
    if (str(alignment_in[i].seq).count('-')/aln_length) < (100 - args.threshold)/100:
      SeqIO.write(alignment_in[i], handle=handle, format="fasta-2line") # append to file if all is good




