library(tidyverse)
library(Biostrings)

repeat_family <- "piggyBac"
query_names <- list.files(paste0("data/Repbase/", repeat_family, "_YW/tbl/"))

for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/Repbase/", repeat_family, "_YW/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/Repbase/", repeat_family, "_YW/tbl/", query_names[i]),
                              col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart", "send",
                                            "evalue", "bitscore", "species", "iteration",
                                            "qlen", "Hit_seq", "Hit_class")))
  } 
}

best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

table(best_hits$Hit_class)

bleeding <- best_hits %>%
  filter(!Hit_class %in% c(repeat_family, "Ginger1", "DNA", "Ginger2/TDD", "MuDR", "EnSpm/CACTA", "Transib", "Kolobok",
                           "Harbinger", "ISL2EU", "Mariner/Tc1", "hAT", "Sola1", "Sola2", "Sola3"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("piggyBac")) 


single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, ":", best_hits$sstart, "-", best_hits$send, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/"))){dir.create(paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/", repeat_family, "_YW.fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- paste0(class_in_q$sseqid, ":", class_in_q$sstart, "-", class_in_q$send, "#", gsub(" ", "_", class_in_q$species))
    if(!dir.exists(paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/others/"))){dir.create(paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/others/"))}
    writeXStringSet(class_in_q_seq, paste0("data/Repbase/", repeat_family, "_YW/extracted_fastas/others/", bleeding_classes$Hit_class[i], "_", repeat_family, "_YW.fasta"))
  }

}
