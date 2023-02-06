library(tidyverse)

taxonomy_info <- read_tsv("data/taxid_probe_2/taxonomy.txt",
                          col_names = c("TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                          show_col_types = F) %>%
  filter(domain %in% c("Bacteria", "Eukaryota", "Archaea", "Viruses")) %>%
  base::unique()

now_found <- read_tsv("data/taxid_probe/compiled_missing_accessions_completed.tsv") %>%
  dplyr::rename(domain = kingdom)

no_phylum <- now_found[now_found$phylum == "-",] %>%
  mutate(kingdom = "-")

kingdom_phylum <- taxonomy_info %>% filter(phylum != "-") %>% dplyr::select(kingdom, phylum) %>% base::unique()

now_found_kingdom_phylum <- now_found %>%
  inner_join(kingdom_phylum) %>%
  rbind(no_phylum) %>%
  dplyr::select(Accession, TaxId, ScientificName, domain, kingdom, phylum, class, order, family, genus) 

write_tsv(now_found_kingdom_phylum, "data/taxid_probe/compiled_missing_accessions_completed_fixed.tsv")
