suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
library(cowplot)

# library(optparse)
# 
# option_list = list(
#   make_option(c("-f", "--family"), type="character", default=NULL, 
#               help="family name (e.g. hAT)", metavar="character")
# )
# 
# opt_parser = OptionParser(option_list=option_list)
# opt = parse_args(opt_parser)

classes <- read_tsv("seq_pass_0/classes_by_size.txt", col_names = "family")

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

resolved_stage_6 <- read_tsv("data/taxid_probe/stage_6/stage_6_found_species.tsv",
                             col_names = c("TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)

resolved_stage_7 <- read_tsv("data/taxid_probe/stage_7/stage_7_found_species.tsv",
                             col_names = c("ScientificName", "TaxId", "NCBIScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F) %>%
  group_by(ScientificName) %>%
  dplyr::slice(1) %>%
  ungroup()

resolved_stage_8 <- read_tsv("data/taxid_probe/stage_8/stage_8_found_species.tsv",
                             col_names = c("TwoWordScientificName", "TaxId", "NCBIScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F) %>%
  group_by(TwoWordScientificName) %>%
  arrange(TaxId) %>%
  dplyr::slice(1) %>%
  ungroup()

resolved_stage_9 <- read_tsv("data/taxid_probe/stage_9/stage_9_found_species.tsv",
                             col_names = c("Accession", "TaxId", "ScientificName", "domain", "kingdom", "phylum", "class", "order", "family", "genus"),
                             show_col_types = F)


# for( i in seq_along(classes$family)){
i=17
opt <- list(family=classes$family[i])

print(opt$family)

# read in and manipulate alignment info
alignment_info <- tibble(centroid = names(Biostrings::readAAStringSet(paste0("data/extracted_fastas/2023_04_threshold/", opt$family, "_realigned_trimmed.fasta")))) %>%
  filter(grepl("#", centroid)) %>%
  dplyr::mutate(Accession = centroid) %>%
  tidyr::separate(Accession, into = c("Accession", "ScientificName"), sep = "#") %>%
  dplyr::mutate(ScientificName = gsub("_", " ", ScientificName))

# read in and manipulate cluster info
cluster_info <- read_tsv(paste0("data/extracted_fastas/cdhit_75/", opt$family, ".fasta.75.cd.clstr.tsv"),
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

  found_info_stage_6 <- still_missing %>%
    dplyr::select(-genus) %>%
    inner_join(resolved_stage_6)
  
  compiled_info <- rbind(compiled_info, found_info_stage_6)
  
}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,] %>% select(Accession, centroid, genus, pident)

if(nrow(still_missing) > 0){
  
  cluster_info_nonsep <- read_tsv(paste0("data/extracted_fastas/cdhit_75/", opt$family, ".fasta.75.cd.clstr.tsv"),
                           col_names = c("ClusterId", "SequenceId", "pident", "length"), show_col_types = F) %>%
    dplyr::mutate(Accession = SequenceId) %>%
    tidyr::separate(Accession, into = c("Accession", "ScientificName", "extra"), sep = "#") %>%
    dplyr::mutate(ScientificName = gsub("_", " ", ScientificName),
                  ScientificName = ifelse(is.na(extra), ScientificName, paste0(ScientificName, "#", extra))) %>%
    dplyr::select(-extra, -pident) %>%
    inner_join(still_missing) %>%
    dplyr::select(-genus)
  
  found_info_stage_7 <- inner_join(cluster_info_nonsep, resolved_stage_7) %>%
    dplyr::select(-NCBIScientificName)
  
  compiled_info <- rbind(compiled_info, found_info_stage_7)

}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,]

if(nrow(still_missing) > 0){

  to_file <- still_missing %>%
    dplyr::select(ScientificName) %>%
    dplyr::mutate(ScientificName = word(ScientificName, 1, 2)) %>%
    base::unique()
  
  write_tsv(to_file, paste0("data/taxid_probe/stage_8/", opt$family, ".txt"), col_names = F)
  
  found_info_stage_8 <- still_missing %>%
    dplyr::mutate(TwoWordScientificName = word(ScientificName, 1, 2)) %>%
    dplyr::select(-genus) %>%
    inner_join(resolved_stage_8) %>%
    dplyr::select(-TwoWordScientificName, -NCBIScientificName)
  
  compiled_info <- rbind(compiled_info, found_info)

}

still_missing <- missing_information[!missing_information$Accession %in% compiled_info$Accession,]

if(nrow(still_missing) > 0){

  found_info_stage_9 <- resolved_stage_9 %>%
    dplyr::select(-genus, -ScientificName) %>%
    inner_join(still_missing)
  
  compiled_info <- rbind(compiled_info, found_info_stage_9)
  
  # write_tsv(tibble(still_missing$Accession), paste0("data/taxid_probe/stage_9/", opt$family, ".txt"), col_names = F)

}

write_tsv(compiled_info, paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv"))

compiled_info <- compiled_info[0,] %>% mutate(superfamily = "")

for(i in 1:18){
  
  in_tbl <- read_tsv(file = paste0("data/taxid_probe/completed_tables/", classes$family[i], ".tsv")) %>% mutate(superfamily = classes$family[i])
  compiled_info <-rbind(compiled_info, in_tbl)
  
}

# compiled_info <- read_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv")) %>%
#   filter(domain == "Eukaryota")
# 
# missing_phyla <- compiled_info %>%
#   filter(phylum != "-") %>%
#   dplyr::select(phylum) %>%
#   base::unique()
# write_tsv(missing_phyla, "data/taxid_probe/protist_phyla.txt", col_names = F)
# 
# missing_class <- compiled_info %>%
#   filter(!phylum %in% missing_phyla$phylum,
#          class != "-") %>%
#   dplyr::select(class) %>%
#   base::unique()
# write_tsv(missing_class, "data/taxid_probe/protist_class.txt", col_names = F)
# 
# missing_order <- compiled_info %>%
#   filter(!phylum %in% missing_phyla$phylum,
#          !class %in% missing_class$class,
#          order != "-") %>%
#   dplyr::select(order) %>%
#   base::unique()
# write_tsv(missing_order, "data/taxid_probe/protist_order.txt", col_names = F)
# 
# missing_family <- compiled_info %>%
#   filter(!phylum %in% missing_phyla$phylum,
#          !class %in% missing_class$class,
#          !order %in% missing_order$order,
#          family != "-") %>%
#   dplyr::select(family) %>%
#   base::unique()
# write_tsv(missing_family, "data/taxid_probe/protist_family.txt", col_names = F)

single_families <- as_tibble(as.data.frame(table(compiled_info$Accession))) %>%
  dplyr::mutate(Var1 = as.character(Var1)) %>%
  dplyr::rename(Accession = Var1) %>%
  dplyr::filter(Freq == 1)

multiple_families <- as_tibble(as.data.frame(table(compiled_info$Accession))) %>%
  dplyr::mutate(Var1 = as.character(Var1)) %>%
  dplyr::rename(Accession = Var1) %>%
  dplyr::filter(Freq >1)

compiled_info <- compiled_info %>%
  mutate(kingdom = ifelse(domain == "Eukaryota" & kingdom == "-", "Protists", kingdom))
  
domain_level <- as_tibble(as.data.frame(table(compiled_info$superfamily, compiled_info$domain))) %>%
  dplyr::filter(as.character(Var2) != "-") %>%
  mutate(superfamily = as.character(Var1), domain  = as.character(Var2))

kingdom_level <- as_tibble(as.data.frame(table(compiled_info[compiled_info$domain == "Eukaryota",]$superfamily, compiled_info[compiled_info$domain == "Eukaryota",]$kingdom))) %>%
  mutate(superfamily = as.character(Var1), kingdom  = as.character(Var2))

combined_level <- domain_level[domain_level$domain != "Eukaryota",] %>%
  dplyr::rename(kingdom=domain) %>%
  rbind(kingdom_level)

# plotting domain level
ggplot(domain_level, aes(x = superfamily, y = Freq, fill = domain)) + geom_bar(position="dodge", stat="identity") +
  facet_grid(domain ~ ., switch = "y", scales = "free", space = "free")+
  theme(strip.text.y.left = element_text(angle = 90),
        legend.position="none") +
  scale_y_continuous(breaks = 0:10*20000) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# plotting kingdom level
ggplot(kingdom_level, aes(x = superfamily, y = Freq, fill = kingdom)) + geom_bar(position="dodge", stat="identity") +
  facet_grid(kingdom ~ ., switch = "y", scales = "free", space = "free")+
  theme(strip.text.y.left = element_text(angle = 90),
        legend.position="none") +
  scale_y_continuous(breaks = 0:10*20000) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# plotting combined maths
combined_level <- inner_join(combined_level, tibble(kingdom = c("Viruses", "Archaea", "Bacteria", "Protists", "Viridiplantae", "Fungi", "Metazoa"), n = 1:7)) %>%
  group_by(superfamily) %>%
  dplyr::mutate(Prop = Freq/sum(Freq)) %>%
  ungroup()

supp.labs <- combined_level$kingdom
names(supp.labs) <- combined_level$n

superfamily_totals <- combined_level %>%
  group_by(superfamily) %>%
  dplyr::mutate(superfamily_sum = sum(Freq)) %>%
  dplyr::ungroup() %>%
  dplyr::select(superfamily, superfamily_sum) %>%
  base::unique()

organism_totals <- combined_level %>%
  group_by(kingdom) %>%
  dplyr::mutate(kingdom_sum = sum(Freq)) %>%
  dplyr::ungroup() %>%
  dplyr::select(kingdom, kingdom_sum) %>%
  base::unique()

# organism pie chart
pie_plot <- ggplot(organism_totals, aes(x="", y=kingdom_sum, fill=kingdom))+
  geom_bar(stat = "identity") +
  geom_col(color = "black") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "null"))

# Line plot above major plot
line_plot <- ggplot(superfamily_totals) +
  geom_point(mapping = aes(x=superfamily, y=superfamily_sum), stat = "identity") +
  geom_line(mapping = aes(x=1:18, y=superfamily_sum)) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # panel.border = element_blank(),
        # panel.grid=element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# standard plot
standard_plot <- ggplot(combined_level, aes(x = superfamily, y = Prop, fill = kingdom)) + geom_bar(position="dodge", stat="identity") +
  facet_grid(n ~ ., switch = "y",
             labeller = labeller(n = supp.labs),
             space = "free") +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  # scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(guide = guide_axis(angle = 90))
standard_plot

# log plot
log_plot <- ggplot(combined_level[combined_level$Freq > 0,], aes(x = superfamily, y = Freq, fill = kingdom)) + geom_bar(position="dodge", stat="identity") +
  facet_grid(n ~ ., switch = "y", space = "free",
             labeller = labeller(n = supp.labs)) +
  theme(legend.position="none",
        # strip.text.y.left = element_text(angle = 90),
        strip.text.y.left = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(trans = "log10") +
  scale_x_discrete(guide = guide_axis(angle = 90))

plot_grid(
  pie_plot, line_plot, NULL, standard_plot,
  ncol = 2,
  align = "v",
  axis = "lr",
  rel_widths = c(1, 4),
  rel_heights = c(1, 6)
)

plot_grid(
  pie_plot, line_plot, NULL, log_plot,
  ncol = 2,
  align = "v",
  axis = "lr",
  rel_widths = c(1, 4),
  rel_heights = c(1, 6)
)




compiled_info_Eukaryota <- compiled_info[compiled_info$domain == "Eukaryota",]

table(compiled_info_Eukaryota[compiled_info_Eukaryota$phylum == "-",]$class)
table(compiled_info_Eukaryota[compiled_info_Eukaryota$class == "Aphelidea",]$phylum)
compiled_info_Eukaryota[compiled_info_Eukaryota$phylum == "-",] %>%
  mutate(case_when(class == "" ~ "",
                   .default = "-"))

phylum_heatmap_tbl <- as_tibble(as.data.frame(table(compiled_info_Eukaryota[compiled_info_Eukaryota$phylum!="-",]$phylum,
                                                    compiled_info_Eukaryota[compiled_info_Eukaryota$phylum!="-",]$superfamily))) %>%
  group_by(Var1) %>%
  mutate(Prop = Freq/sum(Freq)) %>%
  ungroup() %>%
  group_by(Var2) %>%
  mutate(HostProp = Freq/sum(Freq)) %>%
  ungroup()

ggplot(phylum_heatmap_tbl) + geom_tile(aes(x = Var2, y = Var1, fill = Freq)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = F)

ggplot(phylum_heatmap_tbl) + geom_tile(aes(x = Var2, y = Var1, fill = Prop)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = F)

ggplot(phylum_heatmap_tbl) + geom_tile(aes(x = Var2, y = Var1, fill = HostProp)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = F)

phylum_heatmap_tbl_non_ginger <- as_tibble(as.data.frame(table(compiled_info_Eukaryota[compiled_info_Eukaryota$phylum!="-" & compiled_info_Eukaryota$superfamily != "Ginger",]$phylum,
                                                    compiled_info_Eukaryota[compiled_info_Eukaryota$phylum!="-" & compiled_info_Eukaryota$superfamily != "Ginger",]$superfamily))) %>%
  group_by(Var1) %>%
  mutate(Prop = Freq/sum(Freq)) %>%
  ungroup

ggplot(phylum_heatmap_tbl_non_ginger) + geom_tile(aes(x = Var2, y = Var1, fill = Freq)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = F)

ggplot(phylum_heatmap_tbl_non_ginger) + geom_tile(aes(x = Var2, y = Var1, fill = Prop)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = F)


# Further analysis for HTT and genes
# i=13
# opt <- list(family=classes$family[i])
# compiled_info <- read_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv")) %>%
#   filter(domain == "Eukaryota")
# 
# # genes utilising order diversity
# tb=table(compiled_info$ClusterId, compiled_info$order)
# tb3 <- as_tibble(tb, .name_repair = "minimal")
# colnames(tb3) <- c("ClusterId", "order", "Freq")
# tb4 <- spread(tb3, ClusterId, Freq)
# tb5 <- tb4 %>% mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
#   mutate_if(.predicate = is.logical, .funs = as.integer)
# tb6 <- as.data.frame(colSums(tb5[,2:ncol(tb5)])) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column()
# colnames(tb6) <- c("ClusterId", "order")
# order_tb<-tb6[tb6$order > 2,] %>%
#   arrange(-order) %>%
#   mutate(ClusterId = as.integer(ClusterId))
# 
# # genes utilising family diversity
# tb=table(compiled_info$ClusterId, compiled_info$family)
# tb3 <- as_tibble(tb, .name_repair = "minimal")
# colnames(tb3) <- c("ClusterId", "family", "Freq")
# tb4 <- spread(tb3, ClusterId, Freq)
# tb5 <- tb4 %>% mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
#   mutate_if(.predicate = is.logical, .funs = as.integer)
# tb6 <- as.data.frame(colSums(tb5[,2:ncol(tb5)])) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column()
# colnames(tb6) <- c("ClusterId", "family")
# family_tb<-tb6[tb6$family > 2,] %>%
#   arrange(-family) %>%
#   mutate(ClusterId = as.integer(ClusterId))
# 
# # Phylum level HTT check
# # HTT utilising phylum diversity
# tb=table(compiled_info$ClusterId, compiled_info$phylum)
# tb3 <- as_tibble(tb, .name_repair = "minimal")
# colnames(tb3) <- c("ClusterId", "phylum", "Freq")
# tb4 <- spread(tb3, ClusterId, Freq)
# tb5 <- tb4 %>% mutate_if(.predicate = is.numeric, .funs = as.logical) %>%
#   mutate_if(.predicate = is.logical, .funs = as.integer)
# tb6 <- as.data.frame(colSums(tb5[,2:ncol(tb5)])) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column()
# colnames(tb6) <- c("ClusterId", "phylum")
# phylum_tb<-tb6[tb6$phylum > 1,] %>%
#   arrange(-phylum) %>%
#   mutate(ClusterId = as.integer(ClusterId))
# 
# htt_tb <- family_tb[family_tb$ClusterId %in% phylum_tb$ClusterId,]
# gene_tb <- family_tb[!family_tb$ClusterId %in% phylum_tb$ClusterId,]
# 
# htt_list <- base::unique(compiled_info[compiled_info$ClusterId %in% htt_tb$ClusterId,]$centroid)
# gene_list <- base::unique(compiled_info[compiled_info$ClusterId %in% gene_tb$ClusterId,]$centroid)
# 
# normal_list <- base::unique(compiled_info[!compiled_info$ClusterId %in% c(gene_list, htt_list),]$centroid)

# compiled_info_cluster_domain <- compiled_info %>%
#   dplyr::select(ClusterId, domain) %>%
# 
# tb=table(compiled_info_cluster_domain$ClusterId, compiled_info_cluster_domain$domain)
# table(compiled_info_cluster_domain$domain)
# # # Dan's non-tidyverse method
# # tb2 = data.frame("Cluster" = dimnames(tb)[[1]], "Bacteria" = tb[,1] ,"Eukaryota" = tb[,2])
# # compiled_info_cluster_domain_long <- as.data.frame(tb2)
# 
# # My tidyverse method
# tb3 <- as_tibble(tb, .name_repair = "minimal")
# colnames(tb3) <- c("ClusterId", "Domain", "Freq")
# 
# tb4 <- spread(tb3, Domain, Freq)
# 
# # check if any domains have clusters which are both greater than zero and and less than total number of sequences in cluster
# multiple_domains <- tb4[as.logical(rowSums((tb4[,2:ncol(tb4)] < rowSums(tb4[,2:ncol(tb4)]) & tb4[,2:ncol(tb4)] > 0))),]
# 
# compiled_info[compiled_info$ClusterId == "6944",]$ScientificName
# 
# opt$family <- "Ginger"
# compiled_info <- read_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv"))
# compiled_info[compiled_info$domain == "-",]
# 
# to_plot <- as_tibble(as.data.frame(table(compiled_info$domain, compiled_info$kingdom))) %>%
#   filter(Freq > 0) %>%
#   mutate(Var1 = as.character(Var1),
#          Var2 = as.character(Var2),
#          Var2 = ifelse(Var2 == "-", Var1, Var2),
#          Var2 = ifelse(Var2 == "Eukaryota", "Other Eukaryota", Var2),
#          Var1 = ifelse(Var1 == "-", "Synthetic constructs", Var1),
#          Var2 = ifelse(Var2 == "-", "Synthetic constructs", Var2))
# 
# ggplot(to_plot) + geom_bar(mapping = aes(fill = Var2, x = Var1, y = Freq), colour="black", position="stack", stat="identity") +
#   theme_bw() +
#   # scale_fill_brewer(palette="Set1") +
#   scale_x_discrete(expand = c(0,0))  +
#   scale_y_continuous(expand = c(0,0))
#   
# 
# table(compiled_info$kingdom)
# table(compiled_info$phylum)
#  
# # # Get unique information for taxonomy analysis
# # unique_tax_info <- compiled_info %>%
# #   dplyr::select(-ClusterId, -SequenceId, -pident, -length, -Accession, -centroid) %>%
# #   base::unique()
# # 
# # kingdom_info <- as_tibble(as.data.frame(table(unique_tax_info$kingdom)))
# # write_tsv(kingdom_info, paste0("data/taxid_probe/", opt$family, "_kingdom.tsv"), col_names = F)
# # phylum_info <- as_tibble(as.data.frame(table(unique_tax_info$phylum)))
# # write_tsv(phylum_info, paste0("data/taxid_probe/", opt$family, "_phylum.tsv"), col_names = F)
# # class_info <- as_tibble(as.data.frame(table(unique_tax_info$class)))
# # write_tsv(class_info, paste0("data/taxid_probe/", opt$family, "_class.tsv"), col_names = F)
# # order_info <- as_tibble(as.data.frame(table(unique_tax_info$order)))
# # write_tsv(order_info, paste0("data/taxid_probe/", opt$family, "_order.tsv"), col_names = F)
# # 
# # write_tsv(tibble(clusters_centroids$Accession), paste0("data/taxid_probe/", opt$family, "_accessions.tsv"), col_names = F)
# # 
# 
# 
# as_tibble(as.data.frame(table(compiled_info$ClusterId))) %>%
#   arrange(-Freq)

compiled_info
