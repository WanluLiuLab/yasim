library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

source("../../R/lib/featureCounts.R")
source("../../R/lib/yasim.R")

all_fq_stats <- get_fq_stats_data("ce11_badread_nanopore2018.fq.stats", "NANOPORE2018") %>%
    dplyr::inner_join(
        get_fq_stats_data("ce11_badread_nanopore2020.fq.stats", "NANOPORE2020"),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_fq_stats_data("ce11_badread_pacbio2016.fq.stats", "PACBIO2016"),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_fq_stats_data("ce11_pbsim2_r94.fq.stats", "R94"),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_fq_stats_data("ce11_pbsim_clr.fq.stats", "CLR"),
        by = "TRANSCRIPT_ID"
    )

all_featureCounts_data <- get_featureCounts_data(
    "ce11_badread_nanopore2018.fq.bam.fc.tsv",
    "NANOPORE2018",
    "TRANSCRIPT_ID"
) %>%
    dplyr::inner_join(
        get_featureCounts_data(
            "ce11_badread_nanopore2020.fq.bam.fc.tsv",
            "NANOPORE2020",
            "TRANSCRIPT_ID"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_featureCounts_data(
            "ce11_badread_pacbio2016.fq.bam.fc.tsv",
            "PACBIO2016",
            "TRANSCRIPT_ID"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_featureCounts_data(
            "ce11_pbsim2_r94.fq.bam.fc.tsv",
            "R94",
            "TRANSCRIPT_ID"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_featureCounts_data(
            "ce11_pbsim_clr.fq.bam.fc.tsv",
            "CLR",
            "TRANSCRIPT_ID"
        ),
        by = "TRANSCRIPT_ID"
    )


all_depth_data <- get_pileup_stats_merged_data(
    "ce11_badread_nanopore2018_trans.fq.bam.stats.d/pileup_stat.merged.tsv",
    "NANOPORE2018"
) %>%
    dplyr::inner_join(
        get_pileup_stats_merged_data(
            "ce11_badread_nanopore2020_trans.fq.bam.stats.d/pileup_stat.merged.tsv",
            "NANOPORE2020"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_pileup_stats_merged_data(
            "ce11_badread_pacbio2016_trans.fq.bam.stats.d/pileup_stat.merged.tsv",
            "PACBIO2016"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_pileup_stats_merged_data(
            "ce11_pbsim2_r94_trans.fq.bam.stats.d/pileup_stat.merged.tsv",
            "R94"
        ),
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::inner_join(
        get_pileup_stats_merged_data(
            "ce11_pbsim_clr_trans.fq.bam.stats.d/pileup_stat.merged.tsv",
            "CLR"
        ),
        by = "TRANSCRIPT_ID"
    )

transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.gtf.transcripts.tsv",
    show_col_types = FALSE
)


all_data <- transcript_stats %>%
    dplyr::inner_join(all_depth_data, by = "TRANSCRIPT_ID") %>%
    dplyr::inner_join(all_featureCounts_data, by = "TRANSCRIPT_ID") %>%
    dplyr::inner_join(all_fq_stats, by = "TRANSCRIPT_ID") %>%
    dplyr::transmute(
        TRANSCRIPT_ID = TRANSCRIPT_ID,
        GENE_ID = GENE_ID,
        NAIVE_LENGTH = NAIVE_LENGTH,
        TRANSCRIBED_LENGTH = TRANSCRIBED_LENGTH,
        EXON_NUMBER = EXON_NUMBER,
        INPUT_DEPTH = YASIM_NANOPORE2018_INPUT_DEPTH,
        TRANSCRIPTOMICS_ALN_NANOPORE2018_AVG_DEPTH = PILEUP_STATS_YASIM_NANOPORE2018_AVG_DEPTH,
        TRANSCRIPTOMICS_ALN_NANOPORE2020_AVG_DEPTH = PILEUP_STATS_YASIM_NANOPORE2020_AVG_DEPTH,
        TRANSCRIPTOMICS_ALN_PACBIO2016_AVG_DEPTH = PILEUP_STATS_YASIM_PACBIO2016_AVG_DEPTH,
        TRANSCRIPTOMICS_ALN_R94_AVG_DEPTH = PILEUP_STATS_YASIM_R94_AVG_DEPTH,
        TRANSCRIPTOMICS_ALN_CLR_AVG_DEPTH = PILEUP_STATS_YASIM_CLR_AVG_DEPTH,
        FC_NANOPORE2018_AVG_DEPTH = FEATURECOUNTS_NANOPORE2018_ACTUAL_N_OF_READS * 1391.652494 / TRANSCRIBED_LENGTH,
        FC_NANOPORE2020_AVG_DEPTH = FEATURECOUNTS_NANOPORE2018_ACTUAL_N_OF_READS * 1424.785730 / TRANSCRIBED_LENGTH,
        FC_PACBIO2016_AVG_DEPTH = FEATURECOUNTS_PACBIO2016_ACTUAL_N_OF_READS * 1436.293968 / TRANSCRIBED_LENGTH,
        FC_R94_AVG_DEPTH = FEATURECOUNTS_R94_ACTUAL_N_OF_READS * 1013.049713 / TRANSCRIBED_LENGTH,
        FC_CLR_AVG_DEPTH = FEATURECOUNTS_CLR_ACTUAL_N_OF_READS * 1966.099342 / TRANSCRIBED_LENGTH,
        YASIM_NANOPORE2018_ACTUAL_DEPTH = YASIM_NANOPORE2018_ACTUAL_N_OF_READS * 1391.652494 / TRANSCRIBED_LENGTH,
        YASIM_NANOPORE2020_ACTUAL_DEPTH = YASIM_NANOPORE2020_ACTUAL_N_OF_READS * 1424.785730 / TRANSCRIBED_LENGTH,
        YASIM_PACBIO2016_ACTUAL_DEPTH = YASIM_PACBIO2016_ACTUAL_N_OF_READS * 1436.293968 / TRANSCRIBED_LENGTH,
        YASIM_R94_ACTUAL_DEPTH = YASIM_R94_ACTUAL_N_OF_READS * 1013.049713 / TRANSCRIBED_LENGTH,
        YASIM_CLR_ACTUAL_DEPTH = YASIM_CLR_ACTUAL_N_OF_READS * 1966.099342 / TRANSCRIBED_LENGTH,
    )

arrow::write_parquet(all_data, "all_gep_data.parquet")

all_data_long <- all_data %>%
    dplyr::select(!(GENE_ID:EXON_NUMBER)) %>%
    tidyr::gather(key = "Condition", value = "Depth", -TRANSCRIPT_ID)

g <- ggplot(all_data_long) +
    geom_density_ridges_gradient(aes(
        x = Depth + 1,
        y = Condition
    )) +
    scale_x_continuous(
        name = "Depth, log10 transformed",
        trans = "log10"
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Depth of all runs")

ggsave("gep_depth_all.pdf", g, width = 12, height = 8)


all_conditions <- unique(all_data_long$Condition)

#' Sample Distances and Clustering
condition_distance_matrix <- matrix(
    nrow = length(all_conditions),
    ncol = length(all_conditions)
)
colnames(condition_distance_matrix) <- all_conditions
rownames(condition_distance_matrix) <- all_conditions
for (accession1 in all_conditions) {
    for (accession2 in all_conditions) {
        condition_distance_matrix[accession1, accession2] <- dist(
            rbind(
                all_data[[accession1]],
                all_data[[accession2]]
            ),
            method = "euclidean"
        )
    }
}

pheatmap(
    condition_distance_matrix,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename = "gep_condition_level_distance.pdf",
    width = 16,
    height = 12,
    main = "Condition-level GEP difference (Euclidean)"
)
