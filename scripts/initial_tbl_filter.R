suppressMessages(library(tidyverse))
library(optparse)

option_list = list(
  make_option(c("-f", "--family"), type="character", default=NULL, 
              help="family name (e.g. hAT)", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="Directory (e.g. data/nr_pass_1/)", metavar="character"),
  make_option(c("-t", "--tsv"), type="character", default=NULL, 
              help="Table file name (e.g. Mariner-21_CGi#Crassostrea_gigas.tsv)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in blast data
table_1 <- read_tsv(paste0(opt$directory, "/", opt$family, "/tbl/", opt$tsv),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore", "species", "iteration",
                                  "qlen", "stitle", "Hit_seq"),
                    show_col_types = F) %>%
  filter(length/qlen > 0.95, evalue < 1e-5, iteration <=3) %>%
  group_by(sseqid) %>%
  arrange(sseqid, -bitscore) %>%
  dplyr::slice(1) %>%
  ungroup()

# write trimmed table
if(nrow(table_1 > 0)){
  write_tsv(table_1, paste0(opt$directory, "/", opt$family, "/tbl/trimmed/", opt$tsv), col_names = F)
}

