library(tidyverse)
library(Biostrings)

repeat_family <- "Novosib"
first_d = 17
last_e = 35
query_names <- list.files(paste0("data/rearranged/", repeat_family, "_YW/tbl/"))

for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/rearranged/", repeat_family, "_YW/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "stitle", "Hit_seq"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/rearranged/", repeat_family, "_YW/tbl/", query_names[i]),
                              col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart", "send",
                                            "evalue", "bitscore", "species", "iteration",
                                            "qlen", "stitle", "Hit_seq")))
  } 
}


single_species_tbl <- table_1 %>%
  filter(qstart <= first_d, qend >= qlen - last_e) %>%
  group_by(sseqid) %>%
  arrange(sseqid, -length) %>%
  dplyr::slice(1) %>%
  ungroup()

single_species_seq <- AAStringSet(single_species_tbl$Hit_seq)
names(single_species_seq) <- paste0(single_species_tbl$sseqid, ":", single_species_tbl$sstart, "-", single_species_tbl$send, "#", gsub(" ", "_", single_species_tbl$species))
if(!dir.exists(paste0("data/rearranged/", repeat_family, "_YW/extracted_fastas/"))){dir.create(paste0("data/rearranged/", repeat_family, "_YW/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/rearranged/", repeat_family, "_YW/extracted_fastas/", repeat_family, "_YW.fasta"))
