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
parser.add_argument('-g', '--gap_proportion', type=float, default=0.1,
                    help='Maximum proportion of sequence ends allowed to be gaps')
parser.add_argument('-t', '--total_proportion', type=float, default=0.5,
                    help='Maximum proportion of sequences allowed to be gaps')
parser.add_argument('-c', '--conserved_proportion', type=float, default=0.4,
                    help='Minimum proportion of highly conserved aa allowed')
args = parser.parse_args()

# define functions
# end gap counter function
def end_gap_counter(seq):
  end_count=0
  for i, aa in enumerate(seq):
    if(aa == '-'):
      end_count+=1
    else:
      break
  return(end_count)

# read in alignment and determine no of sequences
aln_in=AlignIO.read(args.in_alignment, 'fasta')
aln_len=len(aln_in)

# determine end prop based on args and seq length
end_prop=len(aln_in[0].seq)*args.gap_proportion

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

# write alignment to file
with open(args.out_alignment, 'w') as handle:
  for line in range(len(aln_in)):
    sequence=str(aln_in[line].seq)
    # if first/last is gap, count gaps at ends 
    start_gaps=end_gap_counter(sequence)
    end_gaps=end_gap_counter(sequence[::-1])
    # determine proportion of sequence which is gaps
    total_gaps_prop=sequence.count("-")/len(sequence)
    # set good counter to 0
    good=0
    for d in con:
      if sequence[d['k']] == d['aa']:
        good+=1
    # print(line, aln_in[line].name, good/len(con), start_gaps, end_gaps)
    if good/len(con) > args.conserved_proportion and start_gaps < end_prop and end_gaps < end_prop and total_gaps_prop < args.total_proportion:
      # print(aln_in[line].name, line)
      SeqIO.write(aln_in[line], handle=handle, format="fasta-2line") # append to file if all is good
print(str(len(AlignIO.read(args.out_alignment, 'fasta')))+' of '+str(aln_len)+' sequences preserved')
print()

sequence='VVFSDEGKFRFFGIKGCLVRRKPCTALQKEHIVPTVKHGGDGVMMWG----------CMASNDVGKLVPGPDTVDALKL------------------------NEAPRGSRESLQKCVAWRCPDELFWPVDSRDLNFVFH----------'

start_count=0
for i, aa in enumerate(sequence):
  if(aa == '-'):
    start_count+=1
  else:
    break
print(start_count)
