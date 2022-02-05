library(tidyverse)
salmon_quant_sf <- "transcripts_quant/quant.sf"
yasim_tsv <- "hg38.chr1.gene.depth.tsv"

yasm <- read_tsv(yasim_tsv)
salmon_quant <- read_tsv(salmon_quant_sf) %>% dplyr::filter(NumReads>0)
salmon_quant_transcript_without_ensver <- salmon_quant %>%
    dplyr::mutate(Name=gsub("\\..*$","",Name))

all_table <-  dplyr::inner_join(yasm, salmon_quant_transcript_without_ensver, by=c("gene_name" = "Name"))

g <- ggplot(all_table) + stat_summary(aes(x=depth,y=NumReads))
ggsave("yasim_to_salmon_quant.png", plot=g)

