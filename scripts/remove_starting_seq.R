suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))


# parse input variables
option_list = list(
  make_option(c("-f", "--family"), type="character", default=NULL, 
              help="family name (e.g. hAT)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in in_seq
in_seq <- readAAStringSet(paste0("seq_pass_3/unaligned/", opt$family, ".fasta"))

# read in cluster info
cluster_info <- read_tsv(paste0("seq_pass_3/unaligned/", opt$family, ".fasta.cd.clstr.tsv"), col_names = c("Cluster", "Sequence", "pident", "length"), show_col_types = F)

# determine centroids
cluster_info <- cluster_info[cluster_info$pident == "C",] %>%
  mutate(Centroid = Sequence) %>%
  dplyr::select(Cluster, Centroid) %>%
  inner_join(cluster_info, by = "Cluster")

# identify clusters containing starting seq (Yuean and Wessler etc), but centroid is not in starting seq
starting_not_centroid <- cluster_info %>%
  filter(grepl("#", Centroid),
         !grepl("#", Sequence))

# select centroids only of clusters which do not contain starting sequences
seq_to_add_list <- cluster_info[!cluster_info$Cluster %in% starting_not_centroid$Cluster,] %>%
  dplyr::filter(grepl("#", Centroid)) %>%
  dplyr::select(Centroid) %>%
  base::unique()

# select these sequences from the new sequences
seq_to_add <- in_seq[names(in_seq) %in% seq_to_add_list$Centroid]

# write new sequences to files
writeXStringSet(seq_to_add, paste0("seq_pass_3/unaligned/", opt$family, "_to_add.fasta"))