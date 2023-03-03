#!/usr/bin/env python

import os
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict, Counter
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_alignment', type=str, required=True,
                    help='original_alignment')
parser.add_argument('-o', '--out alignment', type=str, required=True,
                    help='original_alignment')
args = parser.parse_args()

# read in alignment and determine no of sequences
aln_in=AlignIO.read(args.in_alignment, 'fasta')
aln_len=len(aln_in)

# make empty dictionary
r = defaultdict(list)

# go through each sequence
for line in range(len(aln_in)):
  sequence=str(aln_in[line].seq)
  # identify aa at each position and append to appropriate position in dict
  for i, aa in enumerate(list(sequence)):
    r[i].append(aa)

# tally aa at each position
count = {k:Counter(v) for k,v in r.items()}
con=[]
# go through tallied dict
for k, v in count.items():
    # for each position determine if >50% of aa are said aa
    for aa, c in v.items():
        prop = c/aln_len
        # if >50 append with dict
        if prop >= 0.5 and aa != '-':
          con.append({"k": k,"aa" :aa})
print(str(len(con))+" conserved amino acids in alignment of length "+str(len(aln_in[0].seq)))

end_prop=int(len(sequence)/10)

# write alignment to file
with open(args.out_alignment, 'w') as handle:
  for line in range(len(aln_in)):
    sequence=str(aln_in[line].seq)
    # if first/last is gap, count gaps at ends 
    if sequence[0] == '-':
      start_gaps=sequence[0:end_prop].count("-")
    else:
      start_gaps=0
    if sequence[0] == '-':
      end_gaps=sequence[(len(sequence)-end_prop):len(sequence)].count("-")
    else:
      end_gaps=0
    # determine proportion of sequence which is gaps
    total_gaps_prop=sequence.count("-")/len(sequence)
    # set good counter to 0
    good=0
    for d in con:
      if sequence[d['k']] == d['aa']:
        good+=1
    # print(line, aln_in[line].name, good/len(con), start_gaps, end_gaps)
    if good/len(con) > 0.4 and start_gaps < end_prop and end_gaps < end_prop and total_gaps_prop < 0.5:
      # print(aln_in[line].name, line)
      SeqIO.write(aln_in[line], handle=handle, format="fasta-2line") # append to file if all is good
print(str(len(AlignIO.read(args.out_alignment, 'fasta')))+' of '+str(aln_len)+' sequences preserved')
print()
