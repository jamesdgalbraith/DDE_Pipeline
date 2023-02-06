suppressMessages(library(tidyverse))
library(optparse)

option_list = list(
  make_option(c("-f", "--family"), type="character", default=NULL, 
              help="family name (e.g. hAT)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

classes <- read_tsv("classes.txt", col_names = "family")

# read in and filter total taxonomy information
taxonomy_info <- read_tsv("data/taxid_probe/stage_0/taxonomy.txt",
                          col_names = c("TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                          show_col_types = F) %>%
  filter(domain %in% c("Bacteria", "Eukaryota", "Archaea", "Viruses")) %>%
  base::unique()

resolved_stage_1 <- read_tsv("data/taxid_probe/resolved_stage_1.tsv",
                             col_names = c("Accession", "TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

resolved_stage_2 <- read_tsv("data/taxid_probe/stage_2/stage_2_found_species.tsv",
                             col_names = c("TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

resolved_stage_3 <- read_tsv("data/taxid_probe/stage_3/stage_3_found_species.tsv",
                             col_names = c("TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

resolved_stage_4 <- read_tsv("data/taxid_probe/stage_4/stage_4_found_species.tsv",
                             col_names = c("ScientificName", "TaxId", "ActualScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

resolved_stage_5 <- read_tsv("data/taxid_probe/stage_5/stage_5_found_species.tsv",
                             col_names = c("Accession", "TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

for( i in seq_along(classes$family)){
opt <- list(family=classes$family[i])

print(opt$family)

# read in and manipulate alignment info
alignment_info <- tibble(centroid = names(Biostrings::readAAStringSet(paste0("data/extracted_fastas/", opt$family, "_aln.fasta.trimmed")))) %>%
  filter(grepl("#", centroid)) %>%
  dplyr::mutate(Accession = centroid) %>%
  tidyr::separate(Accession, into = c("Accession", "ScientificName"), sep = "#") %>%
  dplyr::mutate(ScientificName = gsub("_", " ", ScientificName))

# read in and manipulate cluster info
cluster_info <- read_tsv(paste0("data/extracted_fastas/", opt$family, ".fasta.75.cd.clstr.tsv"),
                         col_names = c("ClusterId", "SequenceId", "pident", "length"), show_col_types = F) %>%
  dplyr::mutate(Accession = SequenceId) %>%
  tidyr::separate(Accession, into = c("Accession", "ScientificName"), sep = "#") %>%
  dplyr::mutate(ScientificName = gsub("_", " ", ScientificName))

# determine centroids
centroids <- cluster_info %>%
  filter(pident == "C") %>%
  dplyr::select(ClusterId, SequenceId) %>%
  dplyr::rename(centroid = SequenceId)

# compile information
clusters_centroids <- cluster_info %>%
  inner_join(centroids) 
compiled_info <- clusters_centroids %>%
  dplyr::mutate(pident = as.double(ifelse(pident == "C", "100.00", pident))) %>%
  inner_join(taxonomy_info)

# determine any missing taxonomy information
missing_information <- cluster_info %>%
  inner_join(centroids) %>%
  dplyr::mutate(pident = as.double(ifelse(pident == "C", "100.00", pident))) %>%
  filter(!ScientificName %in% taxonomy_info$ScientificName) %>%
  dplyr::mutate(genus = sub(" .*", "", ScientificName))

# add data found in stage 1 probe (used accession number)
if(nrow(missing_information[missing_information$Accession %in% resolved_stage_1$Accession,]) > 0){
  
  found_info_stage_1 <- missing_information[missing_information$Accession %in% resolved_stage_1$Accession,] %>%
    dplyr::select(-ScientificName, -genus) %>%
    inner_join(resolved_stage_1)
  compiled_info <- rbind(compiled_info, found_info_stage_1)
  
}

# add data found in stage 2 probe (used truncated species number)
if(nrow(missing_information[missing_information$Accession %in% resolved_stage_1$Accession,]) > 0){

  found_info_stage_2 <- missing_information[!missing_information$Accession %in% resolved_stage_1$Accession,] %>%
    mutate(ScientificName = sub("^(.*? .*?) .*", "\\1", ScientificName)) %>%
    dplyr::select(-genus) %>%
    inner_join(resolved_stage_2)
  compiled_info <- rbind(compiled_info, found_info_stage_2)

}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,] %>%
  mutate(ScientificName = gsub("#", " ", gsub(" ", "_", sub(" ", "#", sub(" ", "#", ScientificName)))))

if(nrow(still_missing) > 0){
  
  found_info_stage_3 <- still_missing %>%
    dplyr::select(-genus) %>%
    inner_join(resolved_stage_3)
  
  compiled_info <- rbind(compiled_info, found_info_stage_3)
  
}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,]

if(nrow(still_missing) > 0){

  found_info_stage_4 <- still_missing %>%
    dplyr::select(-genus) %>%
    inner_join(resolved_stage_4) %>%
    mutate(ScientificName = ActualScientificName) %>%
    dplyr::select(-ActualScientificName)
  
  compiled_info <- rbind(compiled_info, found_info_stage_4)

}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,]

if(nrow(still_missing) > 0){

  found_info_stage_5 <- still_missing %>%
    dplyr::select(-genus, -ScientificName) %>%
    inner_join(resolved_stage_5)
  
  compiled_info <- rbind(compiled_info, found_info_stage_5)

}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,]

if(nrow(still_missing) > 0){

  write_tsv(tibble(unique(still_missing$Accession)), paste0("data/taxid_probe/stage_6/", opt$family, ".txt"), col_names = F)

}

compiled_info %>%
  arrange(ClusterId, pident) %>%
  write_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv"))

}

classes

opt$family <- "Ginger"
compiled_info <- read_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv"))
compiled_info[compiled_info$domain == "-",]

to_plot <- as_tibble(as.data.frame(table(compiled_info$domain, compiled_info$kingdom))) %>%
  filter(Freq > 0) %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Var2 = ifelse(Var2 == "-", Var1, Var2),
         Var2 = ifelse(Var2 == "Eukaryota", "Other Eukaryota", Var2),
         Var1 = ifelse(Var1 == "-", "Synthetic constructs", Var1),
         Var2 = ifelse(Var2 == "-", "Synthetic constructs", Var2))

ggplot(to_plot) + geom_bar(mapping = aes(fill = Var2, x = Var1, y = Freq), colour="black", position="stack", stat="identity") +
  theme_bw() +
  # scale_fill_brewer(palette="Set1") +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_continuous(expand = c(0,0))
  

table(compiled_info$kingdom)
table(compiled_info$phylum)
 
# # Get unique information for taxonomy analysis
# unique_tax_info <- compiled_info %>%
#   dplyr::select(-ClusterId, -SequenceId, -pident, -length, -Accession, -centroid) %>%
#   base::unique()
# 
# kingdom_info <- as_tibble(as.data.frame(table(unique_tax_info$kingdom)))
# write_tsv(kingdom_info, paste0("data/taxid_probe/", opt$family, "_kingdom.tsv"), col_names = F)
# phylum_info <- as_tibble(as.data.frame(table(unique_tax_info$phylum)))
# write_tsv(phylum_info, paste0("data/taxid_probe/", opt$family, "_phylum.tsv"), col_names = F)
# class_info <- as_tibble(as.data.frame(table(unique_tax_info$class)))
# write_tsv(class_info, paste0("data/taxid_probe/", opt$family, "_class.tsv"), col_names = F)
# order_info <- as_tibble(as.data.frame(table(unique_tax_info$order)))
# write_tsv(order_info, paste0("data/taxid_probe/", opt$family, "_order.tsv"), col_names = F)
# 
# write_tsv(tibble(clusters_centroids$Accession), paste0("data/taxid_probe/", opt$family, "_accessions.tsv"), col_names = F)
# 


as_tibble(as.data.frame(table(compiled_info$ClusterId))) %>%
  arrange(-Freq)
