library(tidyverse)

repeat_family <- "CMC"
query_names <- list.files(paste0("data/nr/", repeat_family, "_YW/tbl/"))

for( i in seq_along(query_names)){
  if(i == 1){
    table_1 <- read_tsv(paste0("data/nr/", repeat_family, "_YW/tbl/", query_names[i]),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                      "gapopen", "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "species", "iteration",
                                      "qlen", "stitle", "Hit_seq"))
  } else{
    table_1 <- rbind(table_1,
                     read_tsv(paste0("data/nr/", repeat_family, "_YW/tbl/", query_names[i]),
                              col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart", "send",
                                            "evalue", "bitscore", "species", "iteration",
                                            "qlen", "stitle", "Hit_seq")))
  } 
}

# select best hits
best_hits <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(sseqid, -length) %>%
  dplyr::slice(1) %>%
  ungroup()

# get sequences
single_species_seq <- AAStringSet(best_hits$Hit_seq)
names(single_species_seq) <- paste0(best_hits$sseqid, ":", best_hits$sstart, "-", best_hits$send, "#", gsub(" ", "_", best_hits$species))
if(!dir.exists(paste0("data/nr/", repeat_family, "_YW/extracted_fastas/"))){dir.create(paste0("data/nr/", repeat_family, "_YW/extracted_fastas/"))}
writeXStringSet(single_species_seq, paste0("data/nr/", repeat_family, "_YW/extracted_fastas/", repeat_family, "_YW.fasta"))

# plotting
new_hits<-table_1 %>%
  filter(length >= 0.9*qlen) %>%
  group_by(sseqid) %>%
  arrange(sseqid, iteration) %>%
  dplyr::slice(1) %>%
  ungroup() 

new_hits<-as_tibble(as.data.frame(table(new_hits$iteration))) %>% mutate(Var1 = as.integer(Var1))

new_hits$cum_sum <- cumsum(new_hits$Freq)

ggplot(new_hits, aes(x = Var1, y = cum_sum)) + geom_line() +
  scale_x_continuous(limits = c(1,6), expand = c(0,0),
                     name = "Iteration") +
  scale_y_continuous(limits = c(0,max(new_hits$cum_sum)*1.05), expand = c(0,0),
                     name = ("Cumulative elements found")) +
  theme_bw()

# only use 1
for_heatmap_a <- table_1 %>%
  filter(length >= 0.9*qlen) %>%
  dplyr::select(qseqid, sseqid) %>%
  base::unique()

for_heatmap_b <- for_heatmap_a %>% dplyr::rename(qseqid_b = qseqid)

for_heatmap_joined <- full_join(for_heatmap_a, for_heatmap_b) %>%
  dplyr::select(-sseqid)

for_heatmap <- as_tibble(as.data.frame(table(for_heatmap_joined$qseqid, for_heatmap_joined$qseqid_b))) %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))

ggplot(for_heatmap, aes(x = Var1, y = Var2, fill = Freq)) + geom_tile() + theme_bw() +
  scale_fill_distiller(palette = "YlGnBu")
