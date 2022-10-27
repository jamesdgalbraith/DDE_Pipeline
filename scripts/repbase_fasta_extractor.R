library(tidyverse)
library(Biostrings)

# first parallel -a classes.txt --bar -j 1 cat ./data/pass_3/{}/tbl/*.tsv ">" ./data/pass_3/{}/tbl/compiled_{}.out

repeat_families <- read_tsv("classes.txt", col_names = "repeat_family")
pass <- "pass_3"
for (j in 1:nrow(repeat_families)){
  if(!repeat_family %in% c("Ginger", "PHIS", "CMC", "Tc1marPlm")){
  repeat_family <- repeat_families$repeat_family[j]
  table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/compiled_", repeat_family, ".out"),
                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                    "gapopen", "qstart", "qend", "sstart", "send",
                                    "evalue", "bitscore", "species", "iteration",
                                    "qlen", "Hit_seq", "Hit_class"))
  
    best_hits <- table_1 %>%
      filter(length >= 0.9*qlen) %>%
      group_by(sseqid) %>%
      arrange(-bitscore) %>%
      dplyr::slice(1) %>%
      ungroup()
    
    table(best_hits$Hit_class)
    
    bleeding <- best_hits %>%
      filter(!Hit_class %in% c(repeat_family))
    
    best_hits <- best_hits  %>%
      filter(Hit_class %in% c(repeat_family))
    
    best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)
    
    single_species_seq <- AAStringSet(best_hits$Hit_seq)
    names(single_species_seq) <- substr(paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species)), 1, 49)
    if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
    writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))
    
    if(nrow(bleeding) > 0){
      
      bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
      bleeding_classes
      for(i in seq_along(bleeding_classes$Hit_class)){
        class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
        class_in_q
        class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
        names(class_in_q_seq) <- substr(paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species)), 1, 49)
        if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
        writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
      }
    
    }
  }
}

#### Ginger ####
repeat_family <- "Ginger"
table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/compiled_", repeat_family, ".out"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore", "species", "iteration",
                                  "qlen", "Hit_seq", "Hit_class"))

best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

table(best_hits$Hit_class)

bleeding <- best_hits %>%
  filter(!Hit_class %in% c("Ginger1", "Ginger2/TDD"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Ginger1", "Ginger2/TDD"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- substr(paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species)), 1, 49)
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- names(class_in_q_seq) <- substr(paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species)), 1, 49)
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}


#### PHIS ####
table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/compiled_", repeat_family, ".out"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore", "species", "iteration",
                                  "qlen", "Hit_seq", "Hit_class"))

best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

table(best_hits$Hit_class)

bleeding <- best_hits %>%
  filter(!Hit_class %in% c("Harbinger", "ISL2EU"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Harbinger", "ISL2EU"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- substr(paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species)), 1, 49)
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- names(class_in_q_seq) <- substr(paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species)), 1, 49)
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}

#### CMC ####
repeat_family <- "CMC"
table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/compiled_", repeat_family, ".out"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore", "species", "iteration",
                                  "qlen", "Hit_seq", "Hit_class"))

best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

table(best_hits$Hit_class)

bleeding <- best_hits %>%
  filter(!Hit_class %in% c("EnSpm/CACTA"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("EnSpm/CACTA"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- substr(paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species)), 1, 49)
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- names(class_in_q_seq) <- substr(paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species)), 1, 49)
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}

#### Mariner/Tc1 ####
repeat_family <- "Tc1marPlm"

print(repeat_family)
table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/compiled_", repeat_family, ".out"),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))

best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

table(best_hits$Hit_class)

bleeding <- best_hits %>%
  filter(!Hit_class %in% c("Mariner/Tc1"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Mariner/Tc1"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- substr(paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species)), 1, 49)
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- names(class_in_q_seq) <- substr(paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species)), 1, 49)
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}
