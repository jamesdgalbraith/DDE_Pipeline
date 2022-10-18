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
parser.add_argument('-f', '--family', type=str, required=True,
                    help='Name of family to be fixed')
args = parser.parse_args()

rb_seq = 'unaligned_seq/'+args.family+'_Repbase.fasta'
align_seq = 'unaligned_seq/'+args.family+'_YW.fasta'
out_seq = 'unaligned_seq/'+args.family+'_cleaned_Repbase.fasta'

alignment_in = AlignIO.read(align_seq, 'fasta')

for x in range(len(alignment_in)):
  if x == 0:
    id_list = [alignment_in[x].id]
  else:
    id_list.append(alignment_in[x].id)

with open(rb_seq, 'r') as handle:
  with open(out_seq, 'w') as briefcase:
    for record in SeqIO.parse(handle, "fasta"):
        if record.name not in id_list:
          SeqIO.write(record, briefcase, "fasta-2line")

