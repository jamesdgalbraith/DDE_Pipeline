#!/usr/bin/env python

from re import sub
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_file', type=str, required=True,
                    help='Input .clstr file')
args = parser.parse_args()

out_file=args.in_file+'.tsv'

with open(args.in_file, 'r') as f:
    lines = f.readlines()
    with open(out_file, 'w') as o:
      for index, line in enumerate(lines):
        if (line[0] == ">"):
          cluster=line.split()[1]
        else:
          seqname=line.split()[2]
          seqlen=sub("aa,", "", line.split()[1])
          if line.split()[3] == "at":
            pident=sub("%", "", line.split()[4])
          else:
            pident="C"
          o.write(cluster+'\t'+sub(">", "", (sub("\\.\\.\\.", "", seqname)))+'\t'+pident+'\t'+seqlen+'\n')
