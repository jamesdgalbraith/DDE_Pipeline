library(tidyverse)
library(Biostrings)

repeat_families <- read_tsv("starting_seq/classes.txt", col_names = "repeat_family")
pass <- "pass_2"
for (j in 1:nrow(repeat_families)){
  repeat_family <- repeat_families$repeat_family[j]
  query_names <- list.files(paste0("data/", pass, "/", repeat_family, "/tbl/"))
  
  if(!repeat_family %in% c("Ginger", "PHIS", "CMC", "Tc1marPlm")){
    print(repeat_family)
    for( i in seq_along(query_names)){
      if(i == 1){
        table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
                            col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                          "gapopen", "qstart", "qend", "sstart", "send",
                                          "evalue", "bitscore", "species", "iteration",
                                          "qlen", "Hit_seq", "Hit_class"))
      } else{
        table_1 <- rbind(table_1,
                         read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
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
      filter(!Hit_class %in% c(repeat_family))
    
    best_hits <- best_hits  %>%
      filter(Hit_class %in% c(repeat_family))
    
    best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)
    
    single_species_seq <- AAStringSet(best_hits$Hit_seq)
    names(single_species_seq) <- paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species))
    if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
    writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))
    
    if(nrow(bleeding) > 0){
      
      bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
      bleeding_classes
      for(i in seq_along(bleeding_classes$Hit_class)){
        class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
        class_in_q
        class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
        names(class_in_q_seq) <- paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species))
        if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
        writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
      }
    
    }
  }
}

#### Ginger ####
repeat_family <- "Ginger"
query_names <- list.files(paste0("data/", pass, "/", repeat_family, "/tbl/"))

print(repeat_family)
for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
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
  filter(!Hit_class %in% c("Ginger1", "Ginger2/TDD"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Ginger1", "Ginger2/TDD"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species))
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}


#### PHIS ####
repeat_family <- "PHIS"
query_names <- list.files(paste0("data/", pass, "/", repeat_family, "/tbl/"))

print(repeat_family)
for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
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
  filter(!Hit_class %in% c("Harbinger", "ISL2EU"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Harbinger", "ISL2EU"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species))
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}

#### CMC ####
repeat_family <- "CMC"
query_names <- list.files(paste0("data/", pass, "/", repeat_family, "/tbl/"))

print(repeat_family)
for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
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
  filter(!Hit_class %in% c("EnSpm/CACTA"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("EnSpm/CACTA"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species))
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}

#### Mariner/Tc1 ####
repeat_family <- "Tc1marPlm"
query_names <- list.files(paste0("data/", pass, "/", repeat_family, "/tbl/"))

print(repeat_family)
for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "Hit_seq", "Hit_class"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/", pass, "/", repeat_family, "/tbl/", query_names[i]),
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
  filter(!Hit_class %in% c("Mariner/Tc1"))

best_hits <- best_hits  %>%
  filter(Hit_class %in% c("Mariner/Tc1"))

best_hits$Hit_seq <- gsub("\\*", "X", best_hits$Hit_seq)

single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))){dir.create(paste0("data/", pass, "/", repeat_family, "/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/", pass, "/", repeat_family, "/extracted_fastas/", repeat_family, ".fasta"))

if(nrow(bleeding) > 0){
  
  bleeding_classes <- as_tibble(as.data.frame(table(bleeding$Hit_class))) %>% dplyr::rename(Hit_class = Var1) %>% dplyr::mutate(Hit_class = as.character(Hit_class))
  bleeding_classes
  for(i in seq_along(bleeding_classes$Hit_class)){
    class_in_q <- bleeding[toupper(bleeding$Hit_class) == toupper(bleeding_classes$Hit_class[i]),]
    class_in_q
    class_in_q_seq <- AAStringSet(class_in_q$Hit_seq)
    names(class_in_q_seq) <- paste0(class_in_q$sseqid, "#", gsub(" ", "_", class_in_q$species))
    if(!dir.exists(paste0("data/", pass, "/others/extracted_fastas"))){dir.create(path = paste0("data/", pass, "/others/extracted_fastas/"), recursive = T)}
    writeXStringSet(class_in_q_seq, filepath = paste0("data/", pass, "/others/extracted_fastas/", sub("/.*", "", bleeding_classes$Hit_class[i]), "_", repeat_family, ".fasta"))
  }
  
}
