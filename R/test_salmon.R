library(tidyverse)
library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--salmon_quant_sf", help = ".sf file produced by salmon quant.", type = "character")
p <- add_argument(p, "--yasim_tsv", help = "Gene depth .tsv file produced by tasim.", type = "character")
argv <- parse_args(p)

salmon_quant_tsv <- argv$salmon_quant_tsv
yasim_tsv <- argv$yasim_tsv

yasm <- read_tsv(yasim_tsv)
salmon_quant <- read_tsv(salmon_quant_sf) %>% dplyr::filter(NumReads>0)

all_table <-  dplyr::inner_join(yasm, salmon_quant, by=c("gene_name" = "Name"))

g <- ggplot(all_table) + stat_summary(aes(x=depth,y=NumReads))
ggsave("yasim_to_salmon_quant.png", plot=g)

