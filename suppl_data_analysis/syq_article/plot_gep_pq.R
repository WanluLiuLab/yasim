library("tidyverse")

gep_inside_gene_isoform_level_variation_filename <- "gep_inside_gene_isoform_level_variation.parquet"

if (file.exists(gep_inside_gene_isoform_level_variation_filename)) {
    gep_inside_gene_isoform_level_variation <- arrow::read_parquet(
        gep_inside_gene_isoform_level_variation_filename
    )
} else {
    quit(status = 1)
}
g <- ggplot(
    gep_inside_gene_isoform_level_variation,
    aes(x = abs(log10(var)), y = Condition, fill = Condition)
) +
    geom_violin() +
    scale_x_continuous("Variance (log 10 Fold Change)", limit = c(0, 5)) +
    # scale_fill_manual(values = color_compl)+
    ggtitle("Distribution of Variation of Isoforms Inside all Genes") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
    )
ggsave(
    "gep_gene_isoform_level_variation_pq.pdf",
    g,
    width = 14, height = 6
)

gep_isoform_level_variation_filename <- "gep_isoform_level_variation.parquet"
if (file.exists(gep_isoform_level_variation_filename)) {
    gep_isoform_level_variation <- arrow::read_parquet(gep_isoform_level_variation_filename)
} else {
    quit(status = 1)
}
g <- ggplot(
    gep_isoform_level_variation,
    aes(x = abs(log10(var)), y = Condition, fill = Condition)
) +
    geom_violin() +
    scale_x_continuous("Variance (log 10 Fold Change)", limit = c(0, 5)) +
    # scale_fill_manual(values = color_compl)+
    ggtitle("Distribution of Variation among Highest and Lowest 5% isoforms") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
    )
ggsave(
    "gep_isoform_level_variation_pq.pdf",
    g,
    width = 14, height = 6
)



gep_gene_level_variation_filename <- "gep_gene_level_variation.parquet"
if (file.exists(gep_gene_level_variation_filename)) {
    gep_gene_level_variation <- arrow::read_parquet(gep_gene_level_variation_filename)
} else {
    quit(status = 1)
}

g <- ggplot(
    gep_gene_level_variation,
    aes(x = abs(log10(var)), y = Condition, fill = Condition)
) +
    geom_violin() +
    scale_x_continuous("Variance (log 10 Fold Change)", limit = c(0, 5)) +
    # scale_fill_manual(values = color_compl)+
    ggtitle("Distribution of Variation among Highest and Lowest 5% Genes") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
    )

ggsave(
    "gep_gene_level_variation_pq.pdf",
    g,
    width = 14, height = 6
)
