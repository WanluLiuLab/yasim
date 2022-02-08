library(tidyverse)
stringtie_quant_tsv <- "stringtie_quant.tsv"
yasim_tsv <- "hg38.chr1.gene.depth.tsv"

yasm <- read_tsv(yasim_tsv)
stringtie_quant <- read_tsv(stringtie_quant_tsv, quote = "\'")

all_table <-  dplyr::inner_join(yasm, stringtie_quant, by=c("gene_name" = "reference_id"))

g <- ggplot(all_table) + stat_summary(aes(x=depth,y=cov))
ggsave("yasim_to_stringtie_quant.png", plot=g)

