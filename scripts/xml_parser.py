#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import re
import pandas as pd
import csv
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_xml', type=str, required=True,
                    help='Input xml to be parsed')
parser.add_argument('-o', '--out', type=str, required=True,
                    help='Parent output folder')
parser.add_argument('-s', '--start', type=int,
                    help='Maximum allowed starting position within the query sequence')
parser.add_argument('-e', '--end', type=int,
                    help='Minimum allowed ending position within the query sequence')
args = parser.parse_args()

# if directories don't exist make them
out_tsv_folder = args.out+'/tbl/'
out_fasta_folder = args.out+'/fasta/'
if os.path.exists(out_tsv_folder) is False:
  os.makedirs(out_tsv_folder)
if os.path.exists(out_fasta_folder) is False:
  os.makedirs(out_fasta_folder) 

# determine query length
if args.start is not None and args.end is not None:
  file = open(args.in_xml, "r")
  for line in file:
    if re.search("Iteration_query-len", line):
      qlen=int(re.sub('.*>', '', re.sub('</.*', '', line)))
      break
# detemine query lengths be good
if args.start is not None and args.end is not None:
  if args.start < 1 or args.end < 1 or args.start > qlen or args.end > qlen:
    sys.exit("Coordinates must be between 1 and the query length")
  elif args.start > args.end:
    sys.exit("Start must be less than end")
  else:
    Trim = True
elif args.start is None and args.end is None:
  Trim = False
elif args.start is not None or args.end is not None:
  sys.exit("If you have one position you need both")

# set output file paths
out_tsv = re.sub(".xml", ".tsv", re.sub(".*/", out_tsv_folder, args.in_xml))
out_fasta = re.sub(".xml", ".fasta", re.sub(".*/", out_fasta_folder, args.in_xml))

# read in data and extract root
tree = ET.parse(args.in_xml)
root = tree.getroot()

with open(out_tsv, 'wt') as out_tsv:
  # with open(out_fasta, 'wt', newline='') as o:
  tsv_writer = csv.writer(out_tsv, delimiter='\t')
  # go through root, picking out necessary data and writing to file 
  for child_0 in root:
    for child_1 in child_0:
      for child_2 in child_1:
        if child_2.tag == "Iteration_iter-num":
          iteration = child_2.text
        if child_2.tag == "Iteration_query-ID":
          Iteration_query_ID = child_2.text
        if child_2.tag == "Iteration_query-len":
          Iteration_query_len = child_2.text
        for child_3 in child_2:
          for child_4 in child_3:
            if child_4.tag == "Hit_id":
              Hit_id = child_4.text
              if "|" in Hit_id:
                Hit_id = Hit_id.split('|')[1]
            if child_4.tag == "Hit_def":
              Hit_def = re.sub(' >.*', '', child_4.text)
              Hit_species = re.sub(r']', '', re.sub(r'.*\[', '', Hit_def))
            if child_4.tag == "Hit_accession":
              Hit_accession = child_4.text
            for hsp in child_4:
              for child_5 in hsp:
                if child_5.tag == "Hsp_bit-score":
                  Hsp_bit_score = re.sub(r'\..*', '', child_5.text)
                if child_5.tag == "Hsp_evalue":
                  Hsp_evalue = float(child_5.text)
                  Hsp_evalue = f'{Hsp_evalue:.2e}'
                if child_5.tag == "Hsp_query-from":
                  Hsp_query_from = child_5.text
                if child_5.tag == "Hsp_query-to":
                  Hsp_query_to = child_5.text
                if child_5.tag == "Hsp_hit-from":
                  Hsp_hit_from = child_5.text
                if child_5.tag == "Hsp_hit-to":
                  Hsp_hit_to = child_5.text
                if child_5.tag == "Hsp_identity":
                  Hsp_identity = child_5.text
                if child_5.tag == "Hsp_positive":
                  Hsp_positive = child_5.text
                if child_5.tag == "Hsp_gaps":
                  Hsp_gaps = child_5.text
                if child_5.tag == "Hsp_align-len":
                  Hsp_align_len = child_5.text
                if child_5.tag == "Hsp_qseq":
                  Hsp_qseq = child_5.text
                if child_5.tag == "Hsp_hseq":
                  Hsp_hseq = child_5.text
                if child_5.tag == "Hsp_midline":
                  Hsp_midline = child_5.text
              perc_iden = round(100*(int(Hsp_identity)/int(Hsp_align_len)),3)
              perc_iden = f'{perc_iden:.3f}'
              mismatch = Hsp_midline.count(" ") + Hsp_midline.count("+") - Hsp_hseq.count("-") - Hsp_qseq.count("-")
              tsv_writer.writerow([Iteration_query_ID, Hit_id, perc_iden, Hsp_align_len, mismatch, Hsp_gaps, Hsp_query_from, Hsp_query_to, Hsp_hit_from, Hsp_hit_to, Hsp_evalue, Hsp_bit_score, Hit_species, iteration, Iteration_query_len, Hit_def, re.sub("-", "", Hsp_hseq)])
              # if Trim is True:
              #   if int(Hsp_query_from) <= args.start and int (Hsp_query_to) >= args.end:
              #     o.write(re.sub(" ", "_", (">"+Hit_id+":"+Hsp_hit_from+"-"+Hsp_hit_to+"#"+Hit_species+"#Iteration "+iteration+'\n')))
              #     o.write(re.sub("-", "", Hsp_hseq)+'\n')
              # else:
              #   if int(Hsp_align_len) >= 0.9 * int(Iteration_query_len):
              #     o.write(re.sub(" ", "_", (">"+Hit_id+":"+Hsp_hit_from+"-"+Hsp_hit_to+"#"+Hit_species+"#Iteration "+iteration+'\n')))
              #     o.write(re.sub("-", "", Hsp_hseq)+'\n')
