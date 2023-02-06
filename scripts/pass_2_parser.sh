#!/bin/bash

echo $A
./scripts/realignment_cleaner.py -a starting_seq/${A}.fasta -u data/pass_3/${A}/extracted_fastas/${A}.fasta -o data/pass_3/${A}/extracted_fastas/cleaned_${A}.fasta
cd-hit -i data/pass_3/${A}/extracted_fastas/cleaned_${A}.fasta -o data/pass_3/${A}/extracted_fastas/cleaned_${A}.fasta.cd -c 0.8
mafft --localpair --thread 128 --add data/pass_3/${A}/extracted_fastas/cleaned_${A}.fasta.cd starting_seq/${A}.fasta > seq_pass_3/${A}.fasta