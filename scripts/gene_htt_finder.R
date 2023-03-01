library(tidyverse)
library(reshape2)

compiled_cluster_info <- read_tsv("data/taxid_probe/compiled_95_cluster_taxonomy.tsv") %>%
  dplyr::filter(!domain == "-") %>%
  mutate(Cluster = paste0(superfamily, "#", ClusterId))

{multiple_seqs <- as_tibble(as.data.frame(table(compiled_cluster_info$Cluster))) %>% filter(Freq > 1)

singletons <- compiled_cluster_info[!compiled_cluster_info$Cluster %in% as.character(multiple_seqs$Var1),]
compiled_cluster_info <- compiled_cluster_info[compiled_cluster_info$Cluster %in% as.character(multiple_seqs$Var1),]

# identify cross domain clusters
compiled_cluster_info_multiple_tbl <- as_tibble(as.data.frame(table(compiled_cluster_info$Cluster, compiled_cluster_info$domain))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::rename(Cluster = Var1, domain = Var2) %>%
  dcast(formula = Cluster ~ domain, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

compiled_cluster_info_multiple_tbl <- compiled_cluster_info_multiple_tbl[rowSums(compiled_cluster_info_multiple_tbl[,(2:ncol(compiled_cluster_info_multiple_tbl))]) > 1,]

cross_domain <- compiled_cluster_info[compiled_cluster_info$Cluster %in% compiled_cluster_info_multiple_tbl$Cluster,]

single_domain <- compiled_cluster_info[!compiled_cluster_info$Cluster %in% cross_domain$Cluster,]

# identify cross kingdom clusters
single_domain_multiple_kingdom_tbl <- as_tibble(as.data.frame(table(single_domain$Cluster, single_domain$kingdom))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, kingdom = Var2) %>%
  dcast(formula = Cluster ~ kingdom, value.var = "Freq") %>%
  as_tibble()  %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_domain_multiple_kingdom_tbl <- single_domain_multiple_kingdom_tbl[rowSums(single_domain_multiple_kingdom_tbl[,(2:ncol(single_domain_multiple_kingdom_tbl))]) > 1,]

cross_kingdom <- single_domain[single_domain$Cluster %in% single_domain_multiple_kingdom_tbl$Cluster,]

single_kingdom <- single_domain[!single_domain$Cluster %in% cross_kingdom$Cluster,]

# identify cross phlyum clusters
single_kingdom_multiple_phylum_tbl <- as_tibble(as.data.frame(table(single_kingdom$Cluster, single_kingdom$phylum))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, phylum = Var2) %>%
  dcast(formula = Cluster ~ phylum, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_kingdom_multiple_phylum_tbl

single_kingdom_multiple_phylum_tbl <- single_kingdom_multiple_phylum_tbl[rowSums(single_kingdom_multiple_phylum_tbl[,(2:ncol(single_kingdom_multiple_phylum_tbl))]) > 1,]

cross_phylum <- single_kingdom[single_kingdom$Cluster %in% single_kingdom_multiple_phylum_tbl$Cluster,]

single_phylum <- single_kingdom[!single_kingdom$Cluster %in% cross_phylum$Cluster,]

# identify cross class clusters
single_phylum_multiple_class_tbl <- as_tibble(as.data.frame(table(single_phylum$Cluster, single_phylum$class))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, class = Var2) %>%
  dcast(formula = Cluster ~ class, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_phylum_multiple_class_tbl

single_phylum_multiple_class_tbl <- single_phylum_multiple_class_tbl[rowSums(single_phylum_multiple_class_tbl[,(2:ncol(single_phylum_multiple_class_tbl))]) > 1,]

cross_class <- single_phylum[single_phylum$Cluster %in% single_phylum_multiple_class_tbl$Cluster,]

single_class <- single_phylum[!single_phylum$Cluster %in% cross_class$Cluster,]

# identify cross order clusters
single_class_multiple_order_tbl <- as_tibble(as.data.frame(table(single_class$Cluster, single_class$order))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, order = Var2) %>%
  dcast(formula = Cluster ~ order, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_class_multiple_order_tbl

single_class_multiple_order_tbl <- single_class_multiple_order_tbl[rowSums(single_class_multiple_order_tbl[,(2:ncol(single_class_multiple_order_tbl))]) > 1,]

cross_order <- single_class[single_class$Cluster %in% single_class_multiple_order_tbl$Cluster,]

single_order <- single_class[!single_class$Cluster %in% cross_order$Cluster,]

# identify cross family clusters
single_order_multiple_family_tbl <- as_tibble(as.data.frame(table(single_order$Cluster, single_order$family))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, family = Var2) %>%
  dcast(formula = Cluster ~ family, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_order_multiple_family_tbl

single_order_multiple_family_tbl <- single_order_multiple_family_tbl[rowSums(single_order_multiple_family_tbl[,(2:ncol(single_order_multiple_family_tbl))]) > 1,]

cross_family <- single_order[single_order$Cluster %in% single_order_multiple_family_tbl$Cluster,]

single_family <- single_order[!single_order$Cluster %in% cross_family$Cluster,]

# identify cross genus clusters
single_family_multiple_genus_tbl <- as_tibble(as.data.frame(table(single_family$Cluster, single_family$genus))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, genus = Var2) %>%
  dcast(formula = Cluster ~ genus, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_family_multiple_genus_tbl

single_family_multiple_genus_tbl <- single_family_multiple_genus_tbl[rowSums(single_family_multiple_genus_tbl[,(2:ncol(single_family_multiple_genus_tbl))]) > 1,]

cross_genus <- single_family[single_family$Cluster %in% single_family_multiple_genus_tbl$Cluster,]

single_genus <- single_family[!single_family$Cluster %in% cross_genus$Cluster,]

# identify cross species clusters
single_genus_multiple_ScientificName_tbl <- as_tibble(as.data.frame(table(single_genus$Cluster, single_genus$ScientificName))) %>%
  dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  dplyr::filter(Var2 != "-") %>%
  dplyr::rename(Cluster = Var1, ScientificName = Var2) %>%
  dcast(formula = Cluster ~ ScientificName, value.var = "Freq") %>%
  as_tibble() %>%
  mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
  mutate_if(.predicate = is.logical, .funs = as.double)

single_genus_multiple_ScientificName_tbl

single_genus_multiple_ScientificName_tbl <- single_genus_multiple_ScientificName_tbl[rowSums(single_genus_multiple_ScientificName_tbl[,(2:ncol(single_genus_multiple_ScientificName_tbl))]) > 1,]

cross_species <- single_genus[single_genus$Cluster %in% single_genus_multiple_ScientificName_tbl$Cluster,]

single_species <- single_genus[!single_genus$Cluster %in% cross_species$Cluster,]
  }

cross_domain
cross_kingdom
cross_phylum
cross_class
cross_order
cross_family
cross_genus
cross_species

single_domain
single_kingdom
single_phylum
single_class
single_order
single_family
single_genus
single_species
singletons

