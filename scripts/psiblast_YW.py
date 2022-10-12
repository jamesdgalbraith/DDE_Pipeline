#!/usr/bin/env python

from Bio import AlignIO, SeqIO
import Bio.Align
import os
from re import sub
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed. Must be fasta with .fasta extension')
args = parser.parse_args()

# set variables
subfamily = sub('.fasta', '', sub('.*/', '', args.in_seq))
out_path = "data/rearranged/"+subfamily
alignment_in = AlignIO.read(args.in_seq, 'fasta')

# check files and folders exists
if(os.path.exists(args.in_seq) == False):
  sys.exit('File not found')

if os.path.exists(out_path) is False:
  os.makedirs(out_path)
if os.path.exists(out_path+'/alignments/') is False:
  os.makedirs(out_path+'/alignments/') 
if os.path.exists(out_path+'/pssms/') is False:
  os.makedirs(out_path+'/pssms/') 
if os.path.exists(out_path+'/xml/') is False:
  os.makedirs(out_path+'/xml/') 
if os.path.exists(out_path+'/fastas/') is False:
  os.makedirs(out_path+'/fastas/') 
if os.path.exists(out_path+'/new_pssms/') is False:
  os.makedirs(out_path+'/new_pssms/') 

# func to find paths
def job_list(x, y):
  pssm_in_path = x+'/pssms/'+y+".pssm"
  new_pssm_out_path = x+'/new_pssms/'+y+".pssm"
  xml_out_path = x+'/xml/'+y+".xml"
  to_run = ([pssm_in_path, xml_out_path, new_pssm_out_path],)
  return(to_run)

# func to run psiblast
def psiblast_p(paths):
  print(sub('.pssm', '', sub('.*/', '', paths[0])))
  os.system('psiblast -in_pssm '+paths[0]+' -db /home/james/databases/nr/nr -out '+paths[1]+' -num_threads 128 -num_iterations 3 -outfmt 5 -max_target_seqs 999999 -evalue 1e-5 -out_pssm '+paths[2])

paths = ()
for x in range(len(alignment_in)):
  paths += job_list(out_path, alignment_in[x].id)

def pool_handeler():
  p = Pool(8)
  p.map(psiblast_p, paths)

pool_handeler()

# for x in range(len(paths)):
#   psiblast_p(paths[x])


