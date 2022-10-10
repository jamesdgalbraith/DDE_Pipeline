#!/usr/bin/env python

import xml.etree.ElementTree as ET
import re
import pandas as pd

tree = ET.parse("data/psiblast_out/Academ_YW_nr_3.xml")
root = tree.getroot()

for child in root:
  print(child.tag, child.attrib)
  a = child

for top in root:
  print(child.tag, child.attrib, child.text)
    if child.tag == "BlastOutput_iterations":
      iterations
      
      print(child.tag, child.attrib, child.text)


for child_0 in root:
  for child_1 in child_0:
    for child_2 in child_1:
      if child_2.tag == "Iteration_iter-num":
        iteration = child_2.text
      if child_2.tag == "Iteration_query-ID":
        Iteration_query_ID = child_2.text
        # print("Iteration "+iteration)
      for child_3 in child_2:
        for child_4 in child_3:
          if child_4.tag == "Hit_def":
            Hit_species = re.sub(r']', '', re.sub(r'.*\[', '', child_4.text))
            # print(Hit_species)
          if child_4.tag == "Hit_accession":
            Hit_accession = child_4.text
            # print(Hit_accession)
          if Hit_accession == "GFS64405":
            for hsp in child_4:
              for child_5 in hsp:
                if iteration == "1":
                  if child_5.tag == "Hsp_bit-score":
                    Hsp_bit_score = child_5.text
                  if Hsp_evalue == "Hsp_evalue":
                    Hsp_hit_to = child_5.text
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
                  
                  print(Iteration_query_ID, Hit_accession, Hsp_hit_from, Hsp_hit_to)
