library(tidyverse)
library(pheatmap)

make_plot <- function(df){
    return (
        ggplot(df)+
    geom_point(
        aes(x=YASIM_1_ACTUAL_TPM, y=TPM, color=Software),
        alpha=0.5, size=0.5
    )+
    theme_bw() +
    geom_abline(slope = 1, intercept = 0) +
    xlim(0, 2) +
    ylim(0, 4)
    )
}


dwgsim_ce11 <- read_tsv("ce11.transcripts_dwgsim_all.tsv")
nanopore2018_ce11 <- read_tsv("ce11.transcripts_badread_nanopore2018_all.tsv")


g_all_above_zero <- dwgsim_ce11 %>% filter(
    FEATURECOUNTS_1_ACTUAL_TPM>0,
    CPPTETGS_1_ACTUAL_TPM>0,
    HTSEQ_COUNT_1_ACTUAL_TPM > 0,
    SALMON_1_ACTUAL_TPM > 0
) %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, SALMON_1_ACTUAL_TPM, HTSEQ_COUNT_1_ACTUAL_TPM,CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

ggsave("all_above_zero_dwgsim.tiff", width=16, height=10)

g_sum_above_zero <- dwgsim_ce11 %>% filter(
    sum(
        FEATURECOUNTS_1_ACTUAL_TPM,
        CPPTETGS_1_ACTUAL_TPM,
        HTSEQ_COUNT_1_ACTUAL_TPM,
        SALMON_1_ACTUAL_TPM
    )>0
) %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, SALMON_1_ACTUAL_TPM, HTSEQ_COUNT_1_ACTUAL_TPM,CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

ggsave("sum_above_zero_dwgsim.tiff", width=16, height=10)

g_nf <- dwgsim_ce11 %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, SALMON_1_ACTUAL_TPM, HTSEQ_COUNT_1_ACTUAL_TPM,CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

corr_matrix <- tibble()
soft_names <-  c("YASIM", "FEATURECOUNTS", "SALMON", "HTSEQ_COUNT", "CPPTETGS")
for(soft_name1 in soft_names){
    for(soft_name2 in soft_names){
        corr_matrix <- dplyr::bind_rows(
            corr_matrix,
            tibble(s1=soft_name1, s2=soft_name2, cor=cor(
                dwgsim_ce11[[sprintf("%s_1_ACTUAL_TPM", soft_name1)]],
                dwgsim_ce11[[sprintf("%s_1_ACTUAL_TPM", soft_name2)]]
            ))
        )
    }
}
spreaded_corr_matrix <- tidyr::spread(corr_matrix, "s1", "cor") %>%
    dplyr::arrange(s2)%>%
    dplyr::select(sort(soft_names)) %>%
    as.matrix()
row.names(spreaded_corr_matrix) <- sort(soft_names)
pheatmap(spreaded_corr_matrix, cluster_rows = FALSE, cluster_cols = FALSE)

ggsave("non_filtered_dwgsim.tiff", width=16, height=10)

make_bias_plot <- function(df, a1, a2, gt){
    df %>% dplyr::select(
        !!a1, !!a2, !!gt
    ) %>%
        ggplot() +
        geom_point(
            aes(x=log2(!!rlang::sym(a1)/!!rlang::sym(gt)), y=log2(!!rlang::sym(a2)/!!rlang::sym(gt))),
             alpha=0.2
        ) + xlim(-10, 10) + ylim(-10, 10) + theme_bw()
}
make_bias_plot(
    dwgsim_ce11,
    "FEATURECOUNTS_1_ACTUAL_TPM",
    "CPPTETGS_1_ACTUAL_TPM",
    "YASIM_1_ACTUAL_TPM"
)
make_bias_plot(
    dwgsim_ce11,
    "HTSEQ_COUNT_1_ACTUAL_TPM",
    "CPPTETGS_1_ACTUAL_TPM",
    "YASIM_1_ACTUAL_TPM"
)
make_bias_plot(
    dwgsim_ce11,
    "SALMON_1_ACTUAL_TPM",
    "CPPTETGS_1_ACTUAL_TPM",
    "YASIM_1_ACTUAL_TPM"
)

g_all_above_zero <- nanopore2018_ce11 %>% filter(
    FEATURECOUNTS_1_ACTUAL_TPM>0,
    CPPTETGS_1_ACTUAL_TPM>0
) %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

ggsave("all_above_zero_nanopore2018.tiff", width=16, height=10)

g_sum_above_zero <- nanopore2018_ce11 %>% filter(
    sum(
        FEATURECOUNTS_1_ACTUAL_TPM,
        CPPTETGS_1_ACTUAL_TPM
    )>0
) %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

ggsave("sum_above_zero_nanopore2018.tiff", width=16, height=10)

g_nf <- nanopore2018_ce11 %>%
     dplyr::select(YASIM_1_ACTUAL_TPM, FEATURECOUNTS_1_ACTUAL_TPM, CPPTETGS_1_ACTUAL_TPM ) %>%
     tidyr::gather(key="Software", value = "TPM", -YASIM_1_ACTUAL_TPM) %>%
    make_plot()

ggsave("non_filtered_nanopore2018.tiff", width=16, height=10)


make_bias_plot(
    nanopore2018_ce11,
    "FEATURECOUNTS_1_ACTUAL_TPM",
    "CPPTETGS_1_ACTUAL_TPM",
    "YASIM_1_ACTUAL_TPM"
)
