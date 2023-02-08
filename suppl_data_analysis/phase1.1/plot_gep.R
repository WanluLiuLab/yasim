library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")
library("glue")
library("rlang")
library("ggrepel")
library("Rtsne")


all_fastq_data <- arrow::read_parquet("all_fastq_data.parquet")

transcript_sample_mean_readlength <- all_fastq_data %>%
    dplyr::group_by(TRANSCRIPT_ID, Condition, TRANSCRIBED_LENGTH) %>%
    dplyr::summarise(
        MEAN_READ_LENGTH = mean(LEN)
    ) %>%
    dplyr::ungroup()


all_fq_stats <- NULL

for (fq_stats_fn in Sys.glob("ce11_*.fq.stats")) {
    condition <- fq_stats_fn %>%
        stringr::str_replace(".fq.stats", "") %>%
        stringr::str_replace("ce11_", "")
    this_fq_stats_data <- readr::read_tsv(
        fq_stats_fn,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            INPUT_DEPTH = col_number(),
            SIMULATED_N_OF_READS = col_number()
        ),
        comment = "#"
    )
    if (is.null(all_fq_stats)) {
        all_fq_stats <- this_fq_stats_data %>%
            dplyr::select(TRANSCRIPT_ID, INPUT_DEPTH) %>%
            dplyr::rename(YASIM_INPUT_DEPTH = INPUT_DEPTH)
    }
    this_mean_length_data <- transcript_sample_mean_readlength %>%
        dplyr::filter(Condition == condition) %>%
        dplyr::select(!Condition)
    this_fq_stats_data <- this_fq_stats_data %>%
        dplyr::select(!(INPUT_DEPTH)) %>%
        dplyr::inner_join(this_mean_length_data, by = "TRANSCRIPT_ID") %>%
        dplyr::mutate(ACTUAL_DEPTH = SIMULATED_N_OF_READS * MEAN_READ_LENGTH / TRANSCRIBED_LENGTH) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("YASIM_%s_ACTUAL_DEPTH", condition)) := ACTUAL_DEPTH
        )
    all_fq_stats <- all_fq_stats %>%
        dplyr::inner_join(this_fq_stats_data, by = "TRANSCRIPT_ID")
    rm(condition, this_fq_stats_data, this_mean_length_data, fq_stats_fn)
}

all_fc_data <- NULL

for (fc_data_fn in Sys.glob("ce11_*.fq.bam.fc.tsv")) {
    condition <- fc_data_fn %>%
        stringr::str_replace(".fq.bam.fc.tsv", "") %>%
        stringr::str_replace("ce11_", "")
    this_fc_data <- readr::read_tsv(
        fc_data_fn,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            Chr = col_character(),
            Start = col_character(),
            End = col_character(),
            Strand = col_character(),
            Length = col_number(),
            NumReads = col_double(),
        ),
        col_names = c(
            "TRANSCRIPT_ID",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
            "NumReads"
        ),
        comment = "#"
    ) %>%
        dplyr::select(TRANSCRIPT_ID, NumReads)
    if (is.null(all_fc_data)) {
        all_fc_data <- this_fc_data %>%
            dplyr::select(TRANSCRIPT_ID)
    }
    this_mean_length_data <- transcript_sample_mean_readlength %>%
        dplyr::filter(Condition == condition) %>%
        dplyr::select(!Condition)
    this_fc_data <- this_fc_data %>%
        dplyr::inner_join(this_mean_length_data, by = "TRANSCRIPT_ID") %>%
        dplyr::mutate(ACTUAL_DEPTH = NumReads * MEAN_READ_LENGTH / TRANSCRIBED_LENGTH) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("FC_%s_ACTUAL_DEPTH", condition)) := ACTUAL_DEPTH
        )
    all_fc_data <- all_fc_data %>%
        dplyr::inner_join(this_fc_data, by = "TRANSCRIPT_ID")
    rm(condition, this_fc_data, this_mean_length_data, fc_data_fn)
}


all_depth_data <- NULL

for (depth_data_fn in Sys.glob("ce11_*_trans.fq.bam.depth.mean.tsv")) {
    condition <- depth_data_fn %>%
        stringr::str_replace("_trans.fq.bam.depth.mean.tsv", "") %>%
        stringr::str_replace("ce11_", "")
    this_depth_data <- readr::read_tsv(
        depth_data_fn,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            AVG_DEPTH = col_number()
        ),
        comment = "#"
    )
    if (is.null(all_depth_data)) {
        all_depth_data <- this_depth_data %>%
            dplyr::select(TRANSCRIPT_ID)
    }
    this_depth_data <- this_depth_data %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("TRANS_%s_ACTUAL_DEPTH", condition)) := AVG_DEPTH
        )
    all_depth_data <- all_depth_data %>%
        dplyr::inner_join(this_depth_data, by = "TRANSCRIPT_ID")
    rm(condition, this_depth_data, depth_data_fn)
}

transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.gtf.transcripts.tsv",
    show_col_types = FALSE
)


all_data <- transcript_stats %>%
    dplyr::full_join(all_depth_data, by = "TRANSCRIPT_ID") %>%
    dplyr::full_join(all_fc_data, by = "TRANSCRIPT_ID") %>%
    dplyr::full_join(all_fq_stats, by = "TRANSCRIPT_ID") %>%
    dplyr::filter(!is.na(.$GENE_ID)) %>%
    replace(is.na(.), 0)

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

ggsave("gep_depth_all.pdf", g, width = 12, height = 18)


g <- ggplot(all_data_long) +
    geom_density_ridges_gradient(aes(
        x = Depth + 1,
        y = Condition
    )) +
    scale_x_continuous(
        name = "Depth, log10 transformed",
        trans = "log10",
        limits = c(1.0001, NA)
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Depth of all runs")

ggsave("gep_depth_all_lim.pdf", g, width = 12, height = 18)

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
    filename = "gep_condition_level_distance.pdf",
    width = 12,
    height = 10,
    main = "Condition-level GEP difference (Euclidean)",
    labels_col = replicate(length(all_conditions), "")
)

all_conditions_pbsims <- as.vector(stringr::str_match(all_conditions, ".*pbsim.*"))
all_conditions_pbsims <- all_conditions_pbsims[!is.na(all_conditions_pbsims)]


condition_distance_matrix_pbsims <- matrix(
    nrow = length(all_conditions_pbsims),
    ncol = length(all_conditions_pbsims)
)
colnames(condition_distance_matrix_pbsims) <- all_conditions_pbsims
rownames(condition_distance_matrix_pbsims) <- all_conditions_pbsims
for (accession1 in all_conditions_pbsims) {
    for (accession2 in all_conditions_pbsims) {
        condition_distance_matrix_pbsims[accession1, accession2] <- dist(
            rbind(
                all_data[[accession1]],
                all_data[[accession2]]
            ),
            method = "euclidean"
        )
    }
}

pheatmap(
    condition_distance_matrix_pbsims,
    filename = "gep_condition_level_distance_pbsims.pdf",
    width = 12,
    height = 10,
    main = "Condition-level GEP difference (Euclidean)",
    labels_col = replicate(length(all_conditions), "")
)

tsne <- Rtsne(all_data[6:length(all_data),], perplexity = 5)
tsne_df <- as.data.frame(tsne$Y)
g <- ggplot(tsne_df) +
    geom_point(aes(
        x = V1,
        y = V2
    )) +
    geom_text_repel(aes(
        x = V1,
        y = V2,
        label = colnames(all_data)[6:length(all_data)]
    )) +
    theme_bw()

ggsave("gep_tsne.pdf", g, width = 16, height = 12)

# Get History of PacBio and ONT
# Look into Phread 33 Quality
# Real data read completeness
