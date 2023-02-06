library(tidyverse)
library(treemap)
library(svglite)
library(VennDiagram)

classes <- read_tsv("classes.txt", col_names = "family")

for(i in seq_along(classes$family)){

opt <- list(family = classes$family[i])
print(opt$family)
compiled_info <- read_tsv(paste0("data/taxid_probe/completed_tables/", opt$family, ".tsv")) %>%
  mutate(superfamily = opt$family)

if(i==1){
    all_compiled_info <- compiled_info
  } else {
    all_compiled_info <- rbind(all_compiled_info, compiled_info)
  }


to_plot <- as_tibble(as.data.frame(table(compiled_info$domain, compiled_info$kingdom))) %>%
  filter(Freq > 0) %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Var2 = ifelse(Var2 == "-", Var1, Var2),
         Var2 = ifelse(Var2 == "Eukaryota", "Other Eukaryota", Var2),
         Var1 = ifelse(Var1 == "-", "Synthetic constructs", Var1),
         Var2 = ifelse(Var2 == "-", "Synthetic constructs", Var2))

plotted_standard <- ggplot(to_plot, aes(x = Var2,  y = Freq, fill = Var2))  +  
  geom_col(position = "dodge") +
  facet_grid(~Var1, scales = "free_x", space = "free_x") +
  theme_bw() +
  # scale_fill_brewer(palette="Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(expand = c(0,0))

ggsave(filename = paste0("plots/", opt$family, "_standard.svg"), plot = plotted_standard, device = "svg", width = 297, height = 210, units = "mm")

plotted_log10 <- ggplot(to_plot, aes(x = Var2,  y = Freq, fill = Var2))  +  
  geom_col(position = "dodge") +
  facet_grid(~Var1, scales = "free_x", space = "free_x") +
  theme_bw() +
  # scale_fill_brewer(palette="Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(trans = "log", breaks = 10^(1:ceiling(log10(max(to_plot$Freq)))),
                     limits = c(1, 10^ceiling(log10(max(to_plot$Freq)))),
                     expand = c(0,0))

ggsave(filename = paste0("plots/", opt$family, "_log10.svg"), plot = plotted_log10, device = "svg", width = 297, height = 210, units = "mm")

# summarise data for treemap
summarised_compiled_info <- compiled_info %>%
  filter(domain !="-") %>%
  mutate(kingdom = ifelse(kingdom == "-" & domain == "Eukaryota", paste0("Other"), kingdom)) %>%
  mutate(kingdom = ifelse(kingdom == "-", "", kingdom)) %>%
  mutate(phylum = ifelse(phylum == "-", paste0("Other"), phylum)) %>%
  select(domain, kingdom, phylum) %>%
  mutate(n=1) %>%
  group_by(domain, kingdom, phylum) %>%
  summarise(Freq = sum(n))

svg(file = paste0("plots/", opt$family, "_treemap.svg"), width = 14.85, height = 10.5)
treemap(summarised_compiled_info,
        index = c("domain", "kingdom", "phylum"),
        vSize = "Freq",
        type="index",
        fontface.labels=c(4,2,1),
        align.labels=list(
          c("left", "top"), 
          c("center", "center"),
          c("right", "bottom")
        ),
        title = paste0("Distribution of ", opt$family, " transposons found in NCBI nr protein database"), aspRatio = 1.414)
dev.off()

}

in_multiple <- read_tsv(paste0("data/taxid_probe/completed_tables/Zator.tsv")) %>%
  mutate(superfamily = "") %>%
  slice(0)

for(i in seq_along(classes$family)){
  
  opt <- list(family = classes$family[i])
  print(opt$family)
  other_compiled_info <- all_compiled_info[all_compiled_info$superfamily != opt$family,]
  family_compiled_info <- all_compiled_info[all_compiled_info$superfamily == opt$family,]
  in_multiple_family <- family_compiled_info[family_compiled_info$SequenceId %in% other_compiled_info$SequenceId,]
  if(nrow(in_multiple_family) >0){
    in_multiple <- rbind(in_multiple, in_multiple_family)
  }
}

100*length(base::unique(in_multiple$SequenceId))/length(base::unique(all_compiled_info$SequenceId))

in_multiple[in_multiple$Accession == "CAH1242273.1",]

in_multiple %>%
  dplyr::select(Accession, domain, kingdom, phylum, class, order, family, genus) %>%
  base::unique()

in_multiple %>%
  dplyr::arrange(Accession)

as_tibble(as.data.frame(table(in_multiple$Accession))) %>%
  mutate(Accession = as.character(Var1)) %>%
  dplyr::select(-Var1) %>%
  inner_join(in_multiple) %>%
  arrange(Accession) %>%
  select(-centroid) %>%
  View()

as_tibble(as.data.frame(table(in_multiple$superfamily)))
as_tibble(as.data.frame(table(all_compiled_info$superfamily)))

all_compiled_info <- all_compiled_info[!is.na(all_compiled_info$Accession),]

CMC_list <- all_compiled_info[all_compiled_info$superfamily == "CMC",]$Accession
Transib_list <- all_compiled_info[all_compiled_info$superfamily == "Transib",]$Accession
other_list <- all_compiled_info[!all_compiled_info$superfamily %in% c("CMC", "Transib"),]$Accession
CMC_Transib_list <- CMC_list[CMC_list %in% Transib_list]

venn.diagram(x = list(CMC_list, Transib_list, other_list),
             category.names = c("CMC" , "Transib " , "Others"),
             filename = 'plots/CMC_Transib.png',
             output=TRUE)

in_multiple_non_cmc_transib <- in_multiple[!in_multiple$Accession %in% CMC_Transib_list,] %>%
  arrange(Accession)
in_multiple_non_cmc_transib %>%
  View()


hAT_list <- all_compiled_info[all_compiled_info$superfamily == "hAT",]$Accession
Ginger_list <- all_compiled_info[all_compiled_info$superfamily == "Ginger",]$Accession
other_list <- all_compiled_info[!all_compiled_info$superfamily %in% c("hAT", "Ginger", "CMC"),]$Accession
hAT_Ginger_list <- hAT_list[hAT_list %in% Ginger_list]

venn.diagram(x = list(hAT_list, Ginger_list, CMC_list, other_list),
             category.names = c("hAT" , "Ginger " , "CMC", "Others"),
             filename = 'plots/hAT_Ginger_CMC.png',
             output=TRUE)

in_multiple_non_cmc_transib <- in_multiple %>%
  filter(Accession %in% hAT_Ginger_list$Accession) %>%
  arrange(Accession)