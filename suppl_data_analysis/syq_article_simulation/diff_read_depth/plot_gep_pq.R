library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
suffix <- args[1]

mutate_confitions <- function(df) {
    df <- df %>% 
    dplyr::filter(
        stringr::str_count(.$Condition, "70") != 0
    )
    return(df)
}

gep_inside_gene_isoform_level_variation_filename <- sprintf(
    "gep_inside_gene_isoform_level_variation_%s.parquet", suffix
)
if (file.exists(gep_inside_gene_isoform_level_variation_filename)) {
    gep_inside_gene_isoform_level_variation <- arrow::read_parquet(
        gep_inside_gene_isoform_level_variation_filename
    ) %>%
    mutate_confitions %>%
    dplyr::mutate(desc="Distribution of Variation of Isoforms Inside all Genes")
} else {
    quit(status=1)
}

gep_isoform_level_variation_filename <- sprintf(
    "gep_isoform_level_variation_%s.parquet", suffix
)
if (file.exists(gep_isoform_level_variation_filename)) {
    gep_isoform_level_variation <- arrow::read_parquet(gep_isoform_level_variation_filename) %>%
    mutate_confitions %>%
    dplyr::mutate(desc="Distribution of Variation among Highest and Lowest 100 Isoforms")
} else {
   quit(status=1)
}
gep_gene_level_variation_filename <- sprintf(
    "gep_gene_level_variation_%s.parquet", suffix
)
if (file.exists(gep_gene_level_variation_filename)) {
    gep_gene_level_variation <- arrow::read_parquet(gep_gene_level_variation_filename) %>%
    mutate_confitions %>%
    dplyr::mutate(desc="Distribution of Variation among Highest and Lowest 100 Genes")
} else {
    quit(status=1)
}

df <- Reduce(
    rbind,
    list(
        gep_inside_gene_isoform_level_variation,
        gep_gene_level_variation,
        gep_isoform_level_variation
    )
)

color_compl<-c("#FBF8CC", "#FDE4CF", "#F1C0E8", "#CFBAF0", "#FFCFD2", "#A3C4F3", "#90DBF4")
g <- ggplot(
    df,
    aes(y = abs(log10(var)), x = Condition, fill = Condition)
) +
    geom_violin() +
    scale_y_continuous("Variance (log 10 Fold Change)", limit = c(0, 5)) +
    scale_fill_manual(values = color_compl)+
    ggtitle("Distribution of Variations") +
    facet_wrap(.~desc) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5)
    )

ggsave(
    sprintf("gep_variation_%s.pdf", suffix),
    g, width = 14, height = 6
)
