#!/usr/bin/python

from re import sub, search
import pandas as pd

family='MuDR'
outFile='/home/james/DDE_Pipeline/data/nr_pass_1/'+family+'_sseqid_species_nonunique.tsv'
inFile='/home/james/DDE_Pipeline/data/nr_pass_1/'+family+'/tbl/compiled_' +family+'.out'

with open(outFile, 'w') as out:
    with open(inFile, 'r') as blast_out:
        for line in blast_out:
            splitline = line.split('\t')
            sseqid = line.split()[1]
            stitle = splitline[15]
            if search('\\[', stitle):
              stitle=stitle.split('[')[1].split(' ')[0:2]
            else:
              stitle=stitle.split(' ')[0:2]
            stitle=' '.join(stitle)
            stitle=sub(']', '', stitle)
            out.write('\t'.join([sseqid, stitle]) + '\n')
