library(tidyverse)
library(arrow)
library(rlang)
library(ggridges)
library(pheatmap)

gene_data <- readr::read_tsv("ce11.ncbiRefSeq.gtf.gene.tsv") %>%
    dplyr::select(GENE_ID, MAPPABLE_LENGTH) %>%
    dplyr::rename(GENE_MAPPABLE_LENGTH = MAPPABLE_LENGTH)

all_gep_data <- gene_data %>%
    dplyr::inner_join(arrow::read_parquet("all_gep_data.parquet"), by="GENE_ID")

gep_data_genelevel <- NULL

for (accesion in colnames(all_gep_data)[11:length(all_gep_data)]){
    message(accesion)
    this_gep_data <- all_gep_data %>%
        dplyr::group_by(GENE_ID) %>%
        dplyr::summarise(
            !!rlang::sym(accesion) := sum(!!rlang::sym(accesion)*TRANSCRIBED_LENGTH)/GENE_MAPPABLE_LENGTH
        ) %>%
        dplyr::ungroup() %>%
        dplyr::distinct()
    print(colnames(this_gep_data))
    if (is.null(gep_data_genelevel)){
        gep_data_genelevel <- this_gep_data
    } else{
        gep_data_genelevel <- gep_data_genelevel %>%
            dplyr::inner_join(this_gep_data, by="GENE_ID")
    }
    rm(accesion, this_gep_data)
}

arrow::write_parquet(gep_data_genelevel, "gep_data_genelevel.parquet")
gep_data_genelevel_filtered <- gep_data_genelevel %>%
    dplyr::mutate(
        sum = rowSums(.[, 2:length(gep_data_genelevel)])
    ) %>%
    dplyr::filter(sum >= 10) %>%
    dplyr::select(!(sum))

gep_data_genelevel_normalized <- gep_data_genelevel_filtered
for (i in 2:length(gep_data_genelevel_filtered)) {
    vec <- gep_data_genelevel_normalized[, i]
    gep_data_genelevel_normalized[, i] <- (vec - min(vec)) / (max(vec) - min(vec)) * 5000
    rm(vec)
}
arrow::write_parquet(
    gep_data_genelevel_normalized,
    "gep_data_genelevel_normalized.parquet"
)

gep_data_genelevel_normalized_long <- gep_data_genelevel_normalized %>%
    tidyr::gather(key="RunAccession", value="Abundance", -GENE_ID)

g <- ggplot(gep_data_genelevel_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance + 1,
        y = RunAccession
    )) +
    scale_x_continuous(
        name = "Normalized NTpG, log10 transformed",
        trans = "log10"
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all runs")
ggsave(
    "scaled_run_level_abundance_genelevel.pdf",
    g,
    width = 8,
    height = 5
)


g <- ggplot(gep_data_genelevel_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance + 1,
        y = RunAccession
    )) +
    scale_x_continuous(
        name = "Normalized NTpG, log10 transformed",
        trans = "log10",
        limits = c(1.0001, 5000)
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all runs")
ggsave(
    "scaled_run_level_abundance_genelevel_lim.pdf",
    g,
    width = 8,
    height = 5
)

all_run_accessions <- colnames(gep_data_genelevel_normalized)[2:length(gep_data_genelevel)]

sample_distance_matrix <- matrix(
    nrow = length(all_run_accessions),
    ncol = length(all_run_accessions)
)
colnames(sample_distance_matrix) <- all_run_accessions
rownames(sample_distance_matrix) <- all_run_accessions
for (accession1 in all_run_accessions) {
    for (accession2 in all_run_accessions) {
        sample_distance_matrix[accession1, accession2] <- dist(
            rbind(
                gep_data_genelevel_normalized[[accession1]],
                gep_data_genelevel_normalized[[accession2]]
            ),
            method = "euclidean"
        )
    }
}

pheatmap(
    sample_distance_matrix,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename = "run_level_distance_genelevel.pdf",
    width = 10,
    height = 8,
    main = "Run-level GEP difference (Euclidean)"
)
