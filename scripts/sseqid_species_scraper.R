library(tidyverse)
library(stringr)

family <- "MuDR"

blast_out <- read_tsv(file = paste0("~/DDE_Pipeline/data/nr_pass_1/compiled_", family, ".out"), 
                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                    "gapopen", "qstart", "qend", "sstart", "send",
                                    "evalue", "bitscore", "species", "iteration",
                                    "qlen", "stitle", "Hit_seq"),
                      show_col_types = F) %>%
  dplyr::select(sseqid, stitle) %>%
  base::unique()
gc()

sseqid_species <- blast_out %>%
  mutate(species = sub(".*?\\[", "", stitle),
         species = word(species, start = 1L, end = 2L, sep = fixed(" ")),
         species = sub("\\]", "", species)) %>%
  dplyr::select(sseqid, species)

write_tsv(x = sseqid_species, file = paste0(family, "_sseqid_species.tsv"))