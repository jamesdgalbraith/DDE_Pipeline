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
parser.add_argument('-d', '--database', type=str, required = True,
                    help='Database to be searched')
parser.add_argument('-o', '--out_dir', type=str, required = True,
                    help='Directory for output')
parser.add_argument('-s', '--db_type', type=str, default = 'prot',
                    help='Sequence type of satabase being searched. Default protein (prot)')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of psiblast processes to run at once')
parser.add_argument('-b', '--blast_threads', type=int, default=1,
                    help='Number of threads to use for each psiblast process')
parser.add_argument('-n', '--iterations', type=int, default=1,
                    help='Number of iterations to run')
args = parser.parse_args()

# set variables
subfamily = sub('.fasta', '', sub('.*/', '', args.in_seq))
out_path = args.out_dir+"/"+subfamily
alignment_in = AlignIO.read(args.in_seq, 'fasta')
folders_to_make=[out_path, out_path+'/split/', out_path+'/alignments/', out_path+'/pssms/', out_path+'/xml/', out_path+'/fasta/', out_path+'/new_pssms/', out_path+'/tbl/', out_path+'/initial_blast/']

# check files and folders exists
if os.path.exists(args.in_seq) == False:
  sys.exit('File not found')

# func to check all out folders exist
def folder_check(x):
  if os.path.exists(x) is False:
    os.makedirs(x)

# func to create pssms and run tblastn/psiblast
def psiblast(x):
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
  os.system('psiblast -subject '+out_path+'/split/'+z+'.fasta -in_msa '+out_path+'/alignments/'+z+'.fasta -out_pssm '+out_path+'/pssms/'+z+'.pssm -out '+out_path+'/initial_blast/'+z+'.out -outfmt 6 -comp_based_stats 1')
  # run blast
  if args.db_type == 'prot':
    os.system('psiblast -in_pssm '+out_path+'/pssms/'+z+'.pssm -db '+args.database+' -out '+out_path+'/xml/'+z+'.xml -num_threads '+str(args.blast_threads)+' -num_iterations '+str(args.iterations)+' -outfmt 5 -max_target_seqs 999999 -evalue 1e-5 -out_pssm '+out_path+'/new_pssms/'+z+'.pssm')
  else:
    os.system('tblastn -in_pssm '+out_path+'/pssms/'+z+'.pssm -db '+args.database+' -out '+out_path+'/xml/'+z+'.xml -num_threads '+str(args.blast_threads)+' -outfmt 5 -max_target_seqs 999999 -evalue 1e-5')

# func to parallelise paths
def pool_handeler():
  p = Pool(args.threads)
  p.map(psiblast, list(range(len(alignment_in))))

# check folders exist
for x in folders_to_make:
  folder_check(x)

# make database
if args.db_type == 'prot':
  if os.path.exists(args.database+'.pdb') == False:
    os.system('makeblastdb -in '+args.database+' -dbtype prot -out '+args.database)
else:
  if os.path.exists(args.database+'.ndb') == False:
    os.system('makeblastdb -in '+args.database+' -dbtype nucl -out '+args.database)

pool_handeler()
