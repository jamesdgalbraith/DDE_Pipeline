suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
library(cowplot)

# Import resolved taxonomy + accession data from 70% clustering information 
compiled_info_70 <- read_tsv("data/taxid_probe/compiled_70_clustering.tsv") %>%
  dplyr::select(Accession, domain, kingdom, phylum, class, order, family, genus) %>%
  base::unique()
# Remove accession information from  70% clustering information 
compiled_info_70_tax_only <- compiled_info_70 %>%
  dplyr::select(-Accession) %>%
  base::unique()
# Import name of classes
classes <- read_tsv("seq_pass_0/classes_by_size.txt", col_names = "family")

compiled_genus_missing <- tibble(Genus = character())

# go along 95% clustering data for each class, identifying accessions and genera not found in 70% clustering information 
for(i in seq_along(classes$family)){
opt <- list(family=classes$family[i])
print(opt$family)
print("read in and manipulate cluster info")
# read in and manipulate cluster info
cluster_info <- read_tsv(paste0("data/extracted_fastas/cdhit_95/", opt$family, ".fasta.95.cd.clstr.tsv"),
                         col_names = c("ClusterId", "SequenceId", "pident", "length"), show_col_types = F) %>%
  dplyr::mutate(Accession = SequenceId) %>%
  tidyr::separate(Accession, into = c("Accession", "ScientificName"), sep = "#") %>%
  dplyr::mutate(ScientificName = gsub("_", " ", ScientificName),
                Genus = sub(" .*", "", ScientificName),
                superfamily = opt$family)
print("compile cluster data")
# compile cluster data
if(i==1){
  all_cluster_info <- cluster_info
  compiled_cluster_info <- cluster_info %>%
    dplyr::select(-Genus) %>%
    inner_join(compiled_info_70)
} else {
  all_cluster_info <- rbind(all_cluster_info, cluster_info)
  compiled_cluster_info <- cluster_info %>%
    dplyr::select(-Genus) %>%
    inner_join(compiled_info_70) %>%
    rbind(compiled_cluster_info)
}
print("identify missing accessions")
# identify missing accessions
accession_missing <- cluster_info[!cluster_info$Accession %in% compiled_info_70$Accession,]
print("identify genera present")
# identify genera present
genus_present <- accession_missing[accession_missing$Genus %in% compiled_info_70_tax_only$genus,]
print("identify genera missing")
# identify genera missing
genus_missing <- accession_missing[!accession_missing$Genus %in% compiled_info_70_tax_only$genus,]

# compile missing genera
print("compile missing genera")
if(i==1){
  all_genus_missing <- genus_missing
} else {
  all_genus_missing <- rbind(all_genus_missing, genus_missing)
}

compiled_genus_missing <- genus_missing %>%
  dplyr::select(Genus) %>%
  rbind(compiled_genus_missing)
}

# identify accessions missing from 70% clustering data
print("identify genera present")
all_accession_missing <- all_cluster_info[!all_cluster_info$Accession %in% compiled_info_70$Accession,] %>%
  dplyr::rename(genus = Genus)

# remove redundant missing genera
compiled_genus_missing <- base::unique(compiled_genus_missing)
# # write missing genera to file for searching of NCBI taxonomy
# write_tsv(compiled_genus_missing, "data/taxid_probe/all_compiled/genus_missing.txt", col_names = F)
# # search NCBI taxonomy usinf esearch in terminal

# read in genera 
compiled_genus_found <- read_tsv("data/taxid_probe/all_compiled/genus_found.tsv", col_names = c("genus", "TaxId", "ActualScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "actual_genus"))

# Remove taxonomy information missing genera
holder <- compiled_genus_found %>%
  dplyr::select(-TaxId, -ActualScientificName) %>%
  base::unique() %>%
  arrange(genus) %>%
  filter(actual_genus != "-")

# Select taxonomy information with only one genus identified (some had many!)
compiled_genus_found_accurate <- as_tibble(as.data.frame(table(holder$genus))) %>%
  mutate(genus = as.character(Var1)) %>%
  filter(Freq == 1)
compiled_genus_found_accurate <- holder[holder$genus%in% compiled_genus_found_accurate$genus,]

found_via_genus <- inner_join(all_accession_missing, compiled_genus_found_accurate) %>%
  dplyr::mutate(genus = actual_genus) %>%
  dplyr::select(-actual_genus)

compiled_cluster_info <- rbind(compiled_cluster_info, found_via_genus)

# Identify accessions of sequences still missing genus information and write to file
missing_accession_post_genus <- all_accession_missing[!all_accession_missing$genus %in% compiled_genus_found_accurate$genus,] %>%
  dplyr::select(Accession) %>%
  base::unique()

# write_tsv(missing_accession_post_genus, "data/taxid_probe/all_compiled/accession_missing.txt", col_names = F)
# Search for accession specific data on NCBI via protein
# Read in identified data
found_accession_post_genus <- read_tsv("data/taxid_probe/all_compiled/accession_found.tsv",
                                       col_names = c("Accession", "TaxId", "ActualScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "actual_genus"))

found_via_accession <- inner_join(all_accession_missing, found_accession_post_genus) %>%
  dplyr::select(-genus, -ScientificName, -TaxId) %>%
  dplyr::rename(genus = actual_genus, ScientificName = ActualScientificName)

compiled_cluster_info <- rbind(compiled_cluster_info, found_via_accession)

still_missing_accession_post_genus <- missing_accession_post_genus[!missing_accession_post_genus$Accession %in% compiled_cluster_info$Accession,]

still_missing <- all_cluster_info[! all_cluster_info$Accession %in% compiled_cluster_info$Accession,] %>%
  dplyr::rename(genus = Genus) %>%
  base::unique()

genus_info <- compiled_cluster_info %>%
  dplyr::select(domain, kingdom, phylum, class, order, family, genus) %>%
  base::unique()

found_via_genus_2 <- inner_join(still_missing, genus_info)

compiled_cluster_info <- rbind(compiled_cluster_info, found_via_genus_2)

write_tsv(compiled_cluster_info, "data/taxid_probe/compiled_95_cluster_taxonomy.tsv")

### FIX BIG PROBLEMS
compiled_cluster_info <- read_tsv("data/taxid_probe/compiled_95_cluster_taxonomy.tsv")

compiled_cluster_info <- compiled_cluster_info %>%
  dplyr::mutate(combined = paste0(SequenceId, "###", superfamily))

determining_multiple <- as_tibble(as.data.frame(table(compiled_cluster_info$combined))) %>%
  dplyr::filter(Freq > 1) %>%
  dplyr::mutate(Var1 = as.character(Var1))

single_compiled_cluster_info <- compiled_cluster_info[!compiled_cluster_info$combined %in% determining_multiple$Var1,] %>%
  dplyr::select(-combined)

# compiled_cluster_info[compiled_cluster_info$combined %in% determining_multiple$Var1,] %>%
#   dplyr::group_by(combined) %>%
#   dplyr::arrange(combined, genus) %>%
#   write_tsv("data/taxid_probe/all_compiled/duplicate_genera.tsv")

fixed_multiple_compiled_cluster_info <- read_tsv("data/taxid_probe/all_compiled/duplicate_genera.tsv") %>%
  dplyr::select(-combined)

compiled_cluster_info <- rbind(single_compiled_cluster_info, fixed_multiple_compiled_cluster_info)

write_tsv(compiled_cluster_info, "data/taxid_probe/compiled_95_cluster_taxonomy.tsv")

