#!/bin/bash

while read A; do
	esearch -db taxonomy -query "${A} [ORGN]" </dev/null | efetch -db taxonomy -format native -mode xml | xtract -pattern Taxon -first TaxId ScientificName -group Taxon -DOMA "(-)" -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" -block "*/Taxon" -match "Rank:superkingdom" -DOMA ScientificName -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName -block "*/Taxon" -match "Rank:class" -CLSS ScientificName -block "*/Taxon" -match "Rank:order" -ORDR ScientificName -block "*/Taxon" -match "Rank:family" -FMLY ScientificName -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName -group Taxon -element "&DOMA" "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" >> taxonomy.txt
done<species.txt

