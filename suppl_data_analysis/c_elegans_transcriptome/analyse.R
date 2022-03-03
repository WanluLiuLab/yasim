library(tidyverse)
library(DESeq2)

all_data <- read_tsv("all.tsv")

deseq_matrix <- all_data %>%
    dplyr::select(
        SALMON_1_ACTUAL_N_OF_READS,
        FEATURECOUNTS_1_ACTUAL_N_OF_READS,
        FEATURECOUNTS_2_ACTUAL_N_OF_READS
        )%>%
    dplyr::filter(rowSums(.) >10) %>%
    dplyr::mutate_if(is.numeric, as.integer) %>%
    as.data.frame()
all_rpkm <- all_data %>%
    dplyr::select(SALMON_1_ACTUAL_RPKM, FEATURECOUNTS_1_ACTUAL_RPKM, FEATURECOUNTS_2_ACTUAL_RPKM) %>%
    tidyr::gather(key="key", value = "value")

ggplot(all_rpkm) + 
    geom_histogram(aes(x=log(value), color=key)) +
    theme_bw()