#!/usr/bin/env python

import os
import re
from os.path import exists
import sys
import argparse

#parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input multi-fasta to be split')
parser.add_argument('-o', '--out_dir', type=str, required=True,
                    help='Output directory')                    
args = parser.parse_args()
out_file=args.out_dir+"unaln_"+args.in_seq.split("/")[-1]
# check file/folder exist
if(exists(args.in_seq) == False):
  sys.exit('File not found')
if(exists(args.out_dir) == False):
  os.mkdir(args.out_dir)

# read file
with open(args.in_seq, 'r') as f:
    lines = f.readlines()
    with open(out_file, 'w') as o:
        for index, line in enumerate(lines):
            if (line[0] == ">"):
                if(index > 0):
                    o.write("\n")
                    o.write(re.sub("", "", line))
                    sys.stdout.write(line)
                else:
                    o.write("\n")
                    o.write(re.sub("-", "", line))
                    # sys.stdout.write(line)
            else:
                line = re.sub("-", "", line).strip()
                o.write(line)
    o.close()
f.close()
