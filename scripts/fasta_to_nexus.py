#!/usr/bin/env python

from Bio import SeqIO, AlignIO
import argparse
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input FASTA MSA')
parser.add_argument('-o', '--out_seq', type=str, required=True,
                    help='Output NEXUS')
args = parser.parse_args()

in_seq = "Academ.fasta"
out_seq = "Academ.nexus"
records = AlignIO.read(in_seq, "fasta")

# AlignIO.write(records, out_seq, "nexus")
# 
# SeqIO.write(records, out_seq, "nexus")

records.annotations[]

record
