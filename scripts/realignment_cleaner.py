#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
from re import sub
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alignment', type=str, required=True,
                    help='Existing alignment')
parser.add_argument('-u', '--unaligned', type=str, required=True,
                    help='Fasta with sequences to be removed')
parser.add_argument('-o', '--out', type=str, required=True,
                    help='Out file')
args = parser.parse_args()

alignment_in = AlignIO.read(args.alignment, 'fasta')

for x in range(len(alignment_in)):
  if x == 0:
    id_list = [alignment_in[x].id]
  else:
    id_list.append(alignment_in[x].id)

with open(args.unaligned, 'r') as handle:
  with open(args.out, 'w') as briefcase:
    for record in SeqIO.parse(handle, "fasta"):
        if record.name not in id_list:
          SeqIO.write(record, briefcase, "fasta-2line")

