library(tidyverse)

cluster_tbl <- read_tsv("data/extracted_fastas/", col_names = c("cluster", "seqnames", "pident", "seqlength"))

cluster_tbl <- cluster_tbl %>%
  dplyr::filter(pident == "C") %>%
  dplyr::rename(centroid = seqnames) %>%
  dplyr::select(cluster, centroid) %>%
  dplyr::inner_join(cluster_tbl) %>%
  mutate(pident = as.double(ifelse(pident == "C", "100.00", pident)),
         cluster = as.character(cluster)) %>%
  tidyr::separate(seqnames, into =c ("acc_no", "species"), sep = "#") %>%
  dplyr::mutate(genus = gsub("^([^_]*)_.*$", "\\1", species)) %>%
  dplyr::mutate(species = gsub("^([^_]*_[^_]*)_.*$", "\\1", species))

species_unique_tbl <- cluster_tbl %>%
  dplyr::select(cluster, species, centroid) %>%
  base::unique()

genus_unique_tbl <- cluster_tbl %>%
  dplyr::select(cluster, genus, centroid) %>%
  base::unique()

genus_count_tbl <- as_tibble(as.data.frame(table(genus_unique_tbl$cluster))) %>%
  filter(Freq > 2) %>%
  arrange(-Freq) %>%
  mutate(cluster = as.character(Var1), genus_no = Freq) %>%
  dplyr::select(cluster, genus_no)

genus_count_tbl

View(cluster_tbl[cluster_tbl$cluster %in% genus_count_tbl$cluster,])

write_tsv(tibble(X1 = cluster_tbl[cluster_tbl$cluster %in% genus_count_tbl$cluster,]$acc_no), "Ginger_multigenus.tsv", col_names = F)
  

