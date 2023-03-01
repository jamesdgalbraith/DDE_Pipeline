library(tidyverse)
library(Biostrings)

tibble_table <- function(x){
  return(dplyr::as_tibble(base::as.data.frame(table(x))))
}

tibble_table_2 <- function(x, y){
  return(dplyr::as_tibble(base::as.data.frame(table(x, y))))
}


taxonomy_info <- read_tsv("data/taxid_probe/compiled_taxonomy.tsv")
superfamilies <- read_tsv("seq_pass_0/classes_by_size.txt", col_names = "superfamily")
compiled_protein_names <- tibble(Accession = character(), Description = character())

compiled_cluster_info <- read_tsv("data/extracted_fastas/cdhit_40/Academ.fasta.40.cd.clstr.tsv", col_names = c("ClusterId", "SequenceId", "pident", "length"))[0,] %>% 
  mutate(superfamily = character())

compiled_protein_names <- tibble(Accession = character(), Description = character())

for(i in seq_along(superfamilies$superfamily)){
  
  cluster_info <- read_tsv(paste0("data/extracted_fastas/cdhit_40/", superfamilies$superfamily[i], ".fasta.40.cd.clstr.tsv"), col_names = c("ClusterId", "SequenceId", "pident", "length"))
  cluster_info$superfamily <- superfamilies$superfamily[i]
  compiled_cluster_info <- rbind(compiled_cluster_info, cluster_info)
  compiled_protein_names <- rbind(read_tsv(paste0("data/protein_names/", superfamilies$superfamily[i], ".tsv")), compiled_protein_names)
}

compiled_cluster_info <- compiled_cluster_info %>%
  mutate(Accession = SequenceId,
         Cluster = paste0(superfamily, "#", ClusterId)) %>%
  tidyr::separate(Accession, into = c("Accession", "ScientificName"), sep = "#") %>%
  inner_join(taxonomy_info)