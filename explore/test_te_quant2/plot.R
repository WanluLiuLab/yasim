library(tidyverse)
df <-readr::read_tsv("all.cmp.tsv.xz")

g <- ggplot(df) +
    geom_bar(aes(y=fn, fill=status), position="stack") +
    theme_bw()
ggsave("plot.pdf", g)
df %>% dplyr::filter(fn=="aln/nhmmer_original.tes.json") %>% head()
