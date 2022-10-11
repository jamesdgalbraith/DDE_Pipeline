#!/usr/bin/env python

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
import Bio.Align
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed. Must be fasta with .fasta extension')
args = parser.parse_args()

# set variables
subfamily = re.sub('.fasta', '', re.sub('.*/', '', args.in_seq))
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


a = list(range(0,len(alignment_in)))
# main loop
for x in a:
  print(alignment_in[x].id)
  # set paths
  alignment_out_path = out_path+'/alignments/'+alignment_in[x].id+".fasta"
  pssm_out_path = out_path+'/pssms/'+alignment_in[x].id+".pssm"
  new_pssm_out_path = out_path+'/new_pssms/'+alignment_in[x].id+".pssm"
  xml_out_path = out_path+'/xml/'+alignment_in[x].id+".xml"
  fasta_out_path = out_path+'/fastas/'+alignment_in[x].id+".fasta"
  # get solo seq
  solo_seq = SeqRecord(seq=Seq(re.sub("-", "", str(alignment_in[x].seq))), id=alignment_in[x].id, description="")
  SeqIO.write(solo_seq, fasta_out_path, "fasta")
  # set up list  for looping
  b = a.copy() # determine seqs to add
  b.remove(x) # remove existing seqs
  # construct alignment
  alignment_out = MultipleSeqAlignment(records=[]) # make empty alignment
  alignment_out.append(alignment_in[x]) # add desired alignment
  for y in b:
    alignment_out.append(alignment_in[y])
  AlignIO.write(alignment_out, alignment_out_path, 'fasta')
  os.system('psiblast -subject '+fasta_out_path+' -in_msa '+alignment_out_path+' -out_pssm '+pssm_out_path+" -out temp")
  os.system('psiblast -in_pssm '+pssm_out_path+' -db /home/james/databases/nr/nr -out '+xml_out_path+' -num_threads 128 -num_iterations 10 -outfmt 5 -max_target_seqs 999999 -evalue 1e-5 -out_pssm '+new_pssm_out_path)
