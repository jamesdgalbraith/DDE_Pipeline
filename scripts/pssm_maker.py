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
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed. Must be fasta with .fasta extension')
parser.add_argument('-o', '--out_dir', type=str, required = True,
                    help='Directory for output')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of psiblast processes to run at once')
args = parser.parse_args()

# set variables
subfamily = sub('.fasta', '', sub('.*/', '', args.in_seq))
out_path = args.out_dir+"/"+subfamily
alignment_in = AlignIO.read(args.in_seq, 'fasta')
folders_to_make=[out_path, out_path+'/split/', out_path+'/alignments/', out_path+'/pssms/', out_path+'/initial_blast/']

# check files and folders exists
if os.path.exists(args.in_seq) == False:
  sys.exit('File not found')

# func to check all out folders exist
def folder_check(x):
  if os.path.exists(x) is False:
    os.makedirs(x)

# func to create pssms and run tblastn/psiblast
def pssm_constructor(x):
  z=alignment_in[x].id # determine sequence name
  print(z)
  b = list(range(len(alignment_in))) # determine seqs to add
  b.remove(x) # remove existing seqs
  alignment_out = MultipleSeqAlignment(records=[]) # make empty alignment
  alignment_out.append(alignment_in[x]) # add desired alignment
  for y in b:
    alignment_out.append(alignment_in[y]) # append all other sequences
  AlignIO.write(alignment_out, out_path+'/alignments/'+z+'.fasta', 'fasta') # write alignment to file
  unaligned=SeqRecord(Seq(sub('-', '', str(alignment_in[x].seq))),id=alignment_in[x].id,description='') # get unaligned sequence
  SeqIO.write(unaligned, out_path+'/split/'+z+'.fasta', 'fasta') # write unaligned sequence to file
  # construct pssm
  os.system('psiblast -subject '+out_path+'/split/'+z+'.fasta -in_msa '+out_path+'/alignments/'+z+'.fasta -out_pssm '+out_path+'/pssms/'+z+'.pssm -out '+out_path+'/initial_blast/'+z+'.out -outfmt 6')
  

# func to parallelise paths
def pool_handeler():
  p = Pool(args.threads)
  p.map(pssm_constructor, list(range(len(alignment_in))))

# check folders exist
for x in folders_to_make:
  folder_check(x)

pool_handeler()
