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

compiled_cluster_info_og <- compiled_cluster_info

compiled_cluster_info <- compiled_cluster_info_og[(compiled_cluster_info_og$SequenceId %in% compiled_cluster_info_og[compiled_cluster_info_og$domain == "-" & compiled_cluster_info_og$pident == "C",]$SequenceId) |
                                                    !compiled_cluster_info_og$SequenceId %in% compiled_cluster_info_og[compiled_cluster_info_og$domain == "-",]$SequenceId,]

multiple_seqs <- as_tibble(as.data.frame(table(compiled_cluster_info$Cluster))) %>% filter(Freq > 1)
singletons <- compiled_cluster_info[!compiled_cluster_info$Cluster %in% as.character(multiple_seqs$Var1),]
compiled_cluster_info <- compiled_cluster_info[compiled_cluster_info$Cluster %in% as.character(multiple_seqs$Var1),]

# identify cross domain clusters
compiled_cluster_info_multiple_tbl <- tibble_table_2(compiled_cluster_info$Cluster, compiled_cluster_info$domain) %>%
  dplyr::mutate(x = as.character(x), y = as.character(y)) %>%
  dplyr::rename(Cluster = x, domain = y) %>%
  dcast(formula = Cluster ~ domain, value.var = "Freq") %>%
  as_tibble() %>%
  dplyr::select(-`-`) %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

compiled_cluster_info_multiple_tbl <- compiled_cluster_info_multiple_tbl[rowSums(compiled_cluster_info_multiple_tbl[,(2:ncol(compiled_cluster_info_multiple_tbl))]) > 1,]
cross_domain <- compiled_cluster_info[compiled_cluster_info$Cluster %in% compiled_cluster_info_multiple_tbl$Cluster,]
single_domain <- compiled_cluster_info[!compiled_cluster_info$Cluster %in% cross_domain$Cluster,]

# count single domain clusters
appendicies <- tibble_table_2(compiled_cluster_info_og$Cluster, compiled_cluster_info_og$domain) %>%
  dplyr::mutate(x = as.character(x), y = as.character(y)) %>%
  dplyr::rename(Cluster = x, domain = y) %>%
  dcast(formula = Cluster ~ domain, value.var = "Freq") %>%
  as_tibble() %>%
  dplyr::select(-`-`) %>%
  group_by(Cluster) %>%
  mutate(appendix = paste0("_", (Archaea + Bacteria + Eukaryota + Viruses)),
         appendix_2 = paste0("_", (Archaea + Bacteria + Eukaryota + Viruses))) %>%
  ungroup() %>%
  dplyr::mutate(appendix = ifelse(Archaea > 0, paste0(appendix, "_", Archaea, "_Archaea"), appendix),
                appendix = ifelse(Bacteria > 0, paste0(appendix, "_", Bacteria, "_Bacteria"), appendix),
                appendix = ifelse(Viruses > 0, paste0(appendix, "_", Viruses, "_Viruses"), appendix),
                appendix = ifelse((Viruses > 0 & Bacteria > 0) |
                                    (Archaea > 0 & Bacteria > 0) |
                                    (Archaea > 0 & Viruses > 0), paste0(appendix, "_multiple"), appendix)) %>%
  dplyr::mutate(appendix_2 = ifelse(Archaea > 0, paste0(appendix_2, "_Archaea", "_", Archaea), appendix_2),
                appendix_2 = ifelse(Bacteria > 0, paste0(appendix_2, "_Bacteria", "_", Bacteria), appendix_2),
                appendix_2 = ifelse(Viruses > 0, paste0(appendix_2, "_Viruses", "_", Viruses), appendix_2),
                appendix_2 = ifelse((Viruses > 0 & Bacteria > 0) |
                                    (Archaea > 0 & Bacteria > 0) |
                                    (Archaea > 0 & Viruses > 0), paste0(appendix_2, "_multiple"), appendix_2)) %>%
  dplyr::select(Cluster, appendix, appendix_2)

appendicies <- compiled_cluster_info_og %>%
  dplyr::filter(pident == "C") %>%
  dplyr::select(Cluster, Accession, superfamily) %>%
  inner_join(appendicies) %>%
  dplyr::mutate(appendix = paste0(Accession, appendix),
                appendix_2 = paste0(Accession, appendix_2)) %>%
  dplyr::select(-Cluster)
 
appendicies[grepl("Bacteria", appendicies$appendix),]
# RENAMING OF SEQUENCES FOR CROSS DOMAIN
for(i in seq_along(superfamilies$superfamily)){
  
  replacement_names <- appendicies[appendicies$superfamily ==superfamilies$superfamily[i],] %>%
    dplyr::select(-superfamily)
  
  write_tsv(replacement_names, paste0("data/extracted_fastas/cdhit_40/gap_check/", superfamilies$superfamily[i], "_domain_nos.txt"), col_names = F)
  
}


# define names for cross superfamily
multiple_family_accessions <- tibble_table(compiled_cluster_info_og$Accession) %>%
  dplyr::filter(Freq > 1) %>%
  dplyr::mutate(x = as.character(x)) %>%
  dplyr::rename(Accession = x)

multiple_family_compiled_cluster_info <- compiled_cluster_info_og%>%
  inner_join(multiple_family_accessions) %>%
  dplyr::group_by(Accession) %>%
  dplyr::mutate(superfamilies = base::paste0(superfamily, collapse = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Accession, superfamily) %>%
  dplyr::mutate(superfamilies = sub(superfamily, "", superfamilies))

single_family_compiled_cluster_info <- compiled_cluster_info_og%>%
  dplyr::filter(!Accession %in% multiple_family_accessions$Accession) %>%
  dplyr::mutate(Freq = 1, superfamilies = "_single")

family_compiled_cluster_info <- base::rbind(single_family_compiled_cluster_info, multiple_family_compiled_cluster_info) %>%
  dplyr::arrange(Cluster)

family_compiled_cluster_info

superfamily_appendicies <- tibble_table_2(family_compiled_cluster_info$Cluster, family_compiled_cluster_info$superfamilies) %>%
  dplyr::mutate(x = as.character(x), y = as.character(y)) %>%
  dplyr::rename(Cluster = x, superfamilies = y) %>%
  filter(Freq > 0) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(Cluster) %>%
  dplyr::mutate(superfamilies = ifelse(superfamilies == "_single", sub("single", "", paste0("_", sub("#.*", "", Cluster))), superfamilies)) %>%
  dplyr::mutate(appendix = base::paste0(superfamilies, "_", Freq, collapse = "")) %>%
  dplyr::ungroup() %>%
  dplyr::select(Cluster, appendix) %>%
    dplyr::arrange(Cluster) %>%
  base::unique()

superfamily_appendicies <- tibble_table(compiled_cluster_info_og$Cluster) %>%
  mutate(Cluster = as.character(x)) %>%
  dplyr::select(-x) %>%
  inner_join(compiled_cluster_info_og) %>%
  dplyr::filter(pident == "C") %>%
  inner_join(superfamily_appendicies) %>%
  mutate(appendix = paste0("_", Freq, appendix))



# RENAMING OF SEQUENCES FOR multiple families
for(i in seq_along(superfamilies$superfamily)){
  
  replacement_names <- superfamily_appendicies[superfamily_appendicies$superfamily ==superfamilies$superfamily[i],] %>%
    mutate(appendix = paste0(Accession, appendix)) %>%
    dplyr::select(Accession, appendix)
  
  write_tsv(replacement_names, paste0("data/extracted_fastas/cdhit_40/gap_check/trees/superfamily_trees/", superfamilies$superfamily[i], "_superfamily_nos.txt"), col_names = F)
  
}
















IS <- read_tsv("data/protein_names/Transib.tsv") %>%
  filter(!grepl("uncharacterized", Description),
         !grepl("Uncharacterized", Description),
         !grepl("unnamed", Description),
         !grepl("predicted protein", Description),
         !grepl("hypothetical", Description)) %>%
  filter(!grepl("KRAB", Description))  %>%
  filter(!grepl("SCAN", Description))  %>%
  filter(!grepl("KRBA", Description)) %>%
  filter(grepl("IS.*family", Description))

tibble_table(sub("family.*", "", IS$Description)) %>%
  arrange(-Freq) %>%
  View()

Tc1_protein_data <- read_tsv("data/nr_pass_1/compiled_Tc1marPlm.out", col_names = F) %>%
  dplyr::select(X2, X16)
colnames(Tc1_protein_data) <- c("Accession", "Description")
Tc1_protein_data <- base::unique(Tc1_protein_data)
Tc1_protein_data[Tc1_protein_data$X2 %in% compiled_cluster_info_og$Accession]

# for(i in seq_along(superfamilies$superfamily)){
#   
#   # read in sequence to rename
#   trimmed_alignment_accession_only <- readAAStringSet(paste0("data/extracted_fastas/cdhit_40/cleaned_alignments/", superfamilies$superfamily[i], "_realigned_trimmed.fasta"))
#   names(trimmed_alignment_accession_only) <- sub("#.*", "", names(trimmed_alignment_accession_only))
#   
#   compiled_cluster_info_og %>%
#     
#   
#   
#   
#   writeXStringSet(trimmed_alignment_accession_only, paste0("data/extracted_fastas/cdhit_40/renamed_alignments/", superfamilies$superfamily[i], ".fasta"))
# }
# 
# cluster_counts <- tibble_table(compiled_cluster_info_og[compiled_cluster_info_og$superfamily == superfamilies$superfamily[i],]$ClusterId) %>%
#   mutate(x = as.double(as.character(x))) %>%
#   dplyr::rename(ClusterId = x)
# cluster_counts <- inner_join(compiled_cluster_info_og[compiled_cluster_info_og$superfamily == superfamilies$superfamily[i],], cluster_counts) %>%
#   dplyr::filter(pident == "C") %>%
#   dplyr::select(Accession, Freq) %>%
#   dplyr::mutate(new_name = paste0(Accession, "___", Freq))
# 
# tibble(seqnames = names(trimmed_alignment_accession_only), start = 1, end = width(trimmed_alignment_accession_only))

tibble_table_2(family_compiled_cluster_info$Cluster, family_compiled_cluster_info$superfamilies) %>%
  dplyr::mutate(x = as.character(x), y = as.character(y)) %>%
  dplyr::rename(Cluster = x, superfamilies = y) %>%
  filter(Freq > 1) %>%
  mutate(Cluster = sub("#.*", "", Cluster)) %>%
  dplyr::select(-Freq) %>%
  base::unique() %>%
  dplyr::arrange(Cluster, superfamilies) %>%
  dplyr::filter(!grepl("single", superfamilies)) %>%
  View()
