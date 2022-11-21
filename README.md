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

`mafft --localpair --thread 128 --add data/pass_2/Sola/extracted_fastas/cleaned_Sola.fasta.cd seq_pass_0/Sola.fasta > seq_pass_2/unclean_Sola.fasta`

ensure sequences retain DDE motifs

`./scripts/DDE_check.py -a seq_pass_0/Sola.fasta -i seq_pass_2/unclean_Sola.fasta -o seq_pass_2/Sola.fasta`

Manually trim alignment in Geneious
