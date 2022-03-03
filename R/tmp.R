library(tidyverse)

all_table <- read_tsv("../makefile_test_data/aaaaaaaaa.tsv")

fig_data_rpkm <- all_table %>%
     dplyr::select(YASIM_MAPPED_1_ACTUAL_RPKM, FEATURECOUNTS_1_ACTUAL_RPKM, CPPTETGS_1_ACTUAL_RPKM) %>%
     tidyr::gather(key="key", value = "value", -YASIM_MAPPED_1_ACTUAL_RPKM)

g <- ggplot(fig_data_rpkm)+
    geom_point(
        aes(x=YASIM_MAPPED_1_ACTUAL_RPKM, y=value, color=key),
        alpha=0.5
    ) +
    geom_abline(slope = 1, intercept = 0)
ggsave("../1_rpkm.png", g, height = 10, width = 15)


fig_data_tpm <- all_table %>%
     dplyr::select(YASIM_MAPPED_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, CPPTETGS_1_ACTUAL_TPM) %>%
     tidyr::gather(key="key", value = "value", -YASIM_MAPPED_1_ACTUAL_TPM)

g <- ggplot(fig_data_tpm)+
    geom_point(
        aes(x=YASIM_MAPPED_1_ACTUAL_TPM, y=value, color=key),
        alpha=0.5
    )+
    geom_abline(slope = 1, intercept = 0)
ggsave("../1_tpm.png", g, height = 10, width = 15)

fig_data_n_reads <- all_table %>%
     dplyr::select(YASIM_MAPPED_1_ACTUAL_N_OF_READS, FEATURECOUNTS_1_ACTUAL_N_OF_READS, CPPTETGS_1_ACTUAL_N_OF_READS) %>%
     tidyr::gather(key="key", value = "value", -YASIM_MAPPED_1_ACTUAL_N_OF_READS)

g <- ggplot(fig_data_n_reads)+
    geom_point(
        aes(x=YASIM_MAPPED_1_ACTUAL_N_OF_READS, y=value, color=key),
        alpha=0.5
    )+
    geom_abline(slope = 1, intercept = 0)
ggsave("../1_n_reads.png", g, height = 10, width = 15)

