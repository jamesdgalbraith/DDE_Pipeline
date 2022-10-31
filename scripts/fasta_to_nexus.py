#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input FASTA MSA')
parser.add_argument('-o', '--out_seq', type=str, required=True,
                    help='Output NEXUS')
args = parser.parse_args()

records = SeqIO.parse(args.in_seq, "fasta")
count = SeqIO.write(records, args.out_seq, "nexus")
print("Converted %i records" % count)
