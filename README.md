# DDE_Pipeline

## Example of initial search run, using second pass of Sola elements as an example:
Search database

`./scripts/initial_psiblast.py -i seq_pass_1/Sola.fasta -b 8 -t 128 -n 3 -d ~/Databases/Repbase/RepBase5May2021.fasta -o data/pass_2 -s nucl`

Convert xml to tbl

`parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/pass_2/Sola/ ::: data/pass_2/Sola/xml/*xml`

compile xmls

`cat data/pass_2/Sola/tbl/* > data/pass_2/Sola/tbl/compiled_Sola.out`

manually extract sequences due to name issues

remove sequences in original alignment

`./scripts/realignment_cleaner.py -a seq_pass_0/Sola.fasta -u data/pass_2/Sola/extracted_fastas/Sola.fasta -o data/pass_2/Sola/extracted_fastas/cleaned_Sola.fasta`

cluster cleaned

`cd-hit -i data/pass_2/Sola/extracted_fastas/cleaned_Sola.fasta -o data/pass_2/Sola/extracted_fastas/cleaned_Sola.fasta.cd -c 0.8`

realign to starting alignment

`mafft --localpair --thread 128 --add data/pass_2/Sola/extracted_fastas/cleaned_Sola.fasta.cd seq_pass_0/Sola.fasta > seq_pass_3/unclean_Sola.fasta`

ensure sequences retain DDE motifs

`./scripts/DDE_check.py -a seq_pass_0/Sola.fasta -i seq_pass_3/unclean_Sola.fasta -o seq_pass_3/Sola.fasta`

Manually trim alignment in Geneious



## Example of NCBI nr search run, using Tc1marPlm elements as an example:

Search database

`./scripts/initial_psiblast.py -i seq_pass_1/Tc1marPlm.fasta -b 8 -t 128 -n 3 -d ~/Databases/nr/nr -o data/nr_pass_1 -s prot`

Convert xml to tbl

`parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/nr_pass_1/Tc1marPlm/ ::: data/nr_pass_1/Tc1marPlm/xml/*xml`

Make directory for trimmed sequences

`mkdir data/nr_pass_1/Tc1marPlm/tbl/trimmed/`

Preprocess BLAST tables

`parallel -j 16 --bar -a Tc1marPlm_list.tsv Rscript scripts/initial_tbl_filter.R -f Tc1marPlm -d data/nr_pass_1/ -t {}`

Compile preprocessed BLAST tables

`cat data/nr_pass_1/Tc1marPlm/tbl/trimmed/*tsv > data/nr_pass_1/Tc1marPlm/tbl/compiled_Tc1marPlm.out`

Extract sequences from table

`Rscript scripts/fasta_extractor_1.R -f Tc1marPlm -d nr_pass_1`

Cluster sequences

`cd-hit -i data/nr_pass_1/Tc1marPlm/extracted_fastas/data/nr_pass_1/Tc1marPlm/extracted_fastas/Tc1marPlm.fasta -o data/nr_pass_1/Tc1marPlm/extracted_fastas/data/nr_pass_1/Tc1marPlm/extracted_fastas/Tc1marPlm.fasta.75.cd -c 0.75 -d 100`

Realign to starting alignment

`mafft --localpair --thread 128 --add data/nr_pass_1/Tc1marPlm/extracted_fastas/data/nr_pass_1/Tc1marPlm/extracted_fastas/Tc1marPlm.fasta.75.cd seq_pass_0/Tc1marPlm.fasta > data/nr_pass_1/Tc1marPlm/extracted_fastas/unclean_Tc1marPlm_aln.fasta`

Ensure sequences retain DDE motifs

`./scripts/DDE_check.py -a seq_pass_0/Tc1marPlm.fasta -i data/nr_pass_1/Tc1marPlm/extracted_fastas/unclean_Tc1marPlm_aln.fasta -o data/extracted_fastas/Tc1marPlm.fasta`
