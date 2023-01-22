library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

#' Get dataset metadata
datasets_metadata <- readr::read_csv(
    "datasets.csv",
    show_col_types = FALSE,
    col_types = c(
        Generation = col_character(),
        SequencerManufacturer = col_character(),
        SequencerModel = col_character(),
        ProjectAccession = col_character(),
        RunAccession = col_character()
    ),
    progress = TRUE
) %>%
    dplyr::mutate(
        gep_filename = sprintf(
            "./datasets/%s/%s.transcript.depth.tsv",
            SequencerModel,
            RunAccession
        )
    )

#' Raw sequence run accessions
all_run_accessions <- unlist(datasets_metadata$RunAccession)
all_sequencer <- unlist(unique(datasets_metadata$SequencerModel))

reference_metadata <- readr::read_tsv(
    "./datasets/reference/ce11_trans.stats",
    show_col_types = FALSE,
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        SEQNAME = col_character(),
        START = col_double(),
        END = col_double(),
        STRAND = col_character(),
        ABSOLUTE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        GC = col_double()
    ),
    progress = TRUE
)


if (file.exists("all_gep_data.parquet")) {
    all_gep_data <- arrow::read_parquet("all_gep_data.parquet")
} else {
    all_gep_data <- reference_metadata
    for (i in seq_len(nrow(datasets_metadata))) {
        data_tuple <- datasets_metadata[i,]
        all_gep_data <-
            readr::read_tsv(
                data_tuple$gep_filename,
                show_col_types = FALSE,
                col_types = c(
                    TRANSCRIPT_ID = col_character(),
                    AVG_DEPTH = col_double()
                ),
                progress = TRUE
            ) %>%
                dplyr::inner_join(all_gep_data, by = "TRANSCRIPT_ID") %>%
                dplyr::mutate(!!rlang::sym(data_tuple$RunAccession) := AVG_DEPTH) %>%
                dplyr::select(!AVG_DEPTH)
        rm(data_tuple)
    }
    arrow::write_parquet(all_gep_data, "all_gep_data.parquet")
}

#' Get number of records
n_of_records <- nrow(all_gep_data)

#' Filter data acording to DESEQ2
all_gep_data_filtered <- all_gep_data %>%
    dplyr::mutate(
        sum = rowSums(.[, 10:length(all_gep_data)])
    ) %>%
    dplyr::filter(sum >= 10) %>%
    dplyr::select(!(sum))
n_of_records <- nrow(all_gep_data_filtered)

#' Get Abundance Information
sample_level_abundance <-
    all_gep_data_filtered[, 10:length(all_gep_data)] %>%
        dplyr::summarise_all(mean) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble(rownames = "RunAccession", .name_repair = "minimal") %>%
        dplyr::rename(Abundance = V1) %>%
        dplyr::arrange(desc(Abundance)) %>%
        dplyr::inner_join(
            datasets_metadata,
            by = "RunAccession"
        )
g <- ggplot(sample_level_abundance) +
    geom_bar(aes(
        y = reorder(RunAccession, Abundance),
        x = Abundance,
        fill = SequencerModel
    ), stat = "identity") +
    theme_bw() +
    scale_x_continuous(
        name = "Mean of Unnormalized NTpI Abundance",
        labels = scales::label_number()
    ) +
    ylab("Run Accession") +
    ggtitle("Bar Plot of Run Accession-Sum of Unnormalized NTpI Abundance")
ggsave(
    "run_level_abundance.pdf",
    g,
    width = 8,
    height = 5
)

all_gep_data_per_base <- all_gep_data_filtered
all_gep_data_per_base[, 10:length(all_gep_data)] <- all_gep_data_filtered[, 10:length(all_gep_data)] * all_gep_data_filtered$TRANSCRIBED_LENGTH / sum(all_gep_data_filtered$TRANSCRIBED_LENGTH)

#' Get Per-Base Abundance Information
sample_level_abundance_per_base <-
    all_gep_data_per_base[, 10:length(all_gep_data)] %>%
        dplyr::summarise_all(sum) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble(rownames = "RunAccession", .name_repair = "minimal") %>%
        dplyr::rename(Abundance = V1) %>%
        dplyr::arrange(desc(Abundance)) %>%
        dplyr::inner_join(
            datasets_metadata,
            by = "RunAccession"
        )
g <- ggplot(sample_level_abundance_per_base) +
    geom_bar(aes(
        y = reorder(RunAccession, Abundance),
        x = Abundance,
        fill = SequencerModel
    ), stat = "identity") +
    theme_bw() +
    scale_x_continuous(
        name = "Mean of Unnormalized Coverage",
        labels = scales::label_number()
    ) +
    ylab("Run Accession") +
    ggtitle("Bar Plot of Run Accession-Sum of Unnormalized Mean Coverage")
ggsave(
    "run_level_abundance_per_base.pdf",
    g,
    width = 8,
    height = 5
)

#' Normalize data by scaling them into (1, 5000)
all_gep_data_normalized <- all_gep_data_filtered
for (i in 10:length(all_gep_data_filtered)) {
    vec <- all_gep_data_normalized[, i]
    all_gep_data_normalized[, i] <- (vec - min(vec)) / (max(vec) - min(vec)) * 5000 + 1
    rm(vec)
}
arrow::write_parquet(
    all_gep_data_normalized,
    "all_gep_data_normalized.parquet"
)

#' Transform into long form
all_gep_data_normalized_long <- all_gep_data_normalized %>%
    dplyr::select(!(GENE_ID:GC)) %>%
    tidyr::gather(
        key = "RunAccession",
        value = "Abundance",
        -TRANSCRIPT_ID
    ) %>%
    dplyr::inner_join(
        datasets_metadata,
        by = "RunAccession"
    )
all_gep_data_filtered_long <- all_gep_data_filtered %>%
    dplyr::select(!(GENE_ID:GC)) %>%
    tidyr::gather(
        key = "RunAccession",
        value = "Abundance",
        -TRANSCRIPT_ID
    ) %>%
    dplyr::inner_join(
        datasets_metadata,
        by = "RunAccession"
    )

model_sorted_all_gep_data_normalized_long <- all_gep_data_normalized_long %>%
    dplyr::mutate(SequencerModel = factor(
        SequencerModel,
        level = sort(unique(.$SequencerModel))
    )) %>%
    dplyr::arrange(desc(SequencerModel)) %>%
    dplyr::mutate(RunAccession = factor(
        RunAccession,
        level = unique(.$RunAccession)
    ))

g <- ggplot(model_sorted_all_gep_data_normalized_long) +
    geom_boxplot(
        aes(
            y = RunAccession,
            x = Abundance,
            color = SequencerModel
        ),
        outlier.alpha = 0.1
    ) +
    theme_bw() +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10"
    ) +
    ylab("Run Accession") +
    ggtitle("Box Plot of Run Accession-NTpI after Normalization")
ggsave(
    "scaled_run_level_abundance.pdf",
    g,
    width = 8,
    height = 5
)

g <- ggplot(model_sorted_all_gep_data_normalized_long) +
    geom_boxplot(
        aes(
            y = RunAccession,
            x = Abundance,
            color = SequencerModel
        ),
        outlier.alpha = 0.1
    ) +
    theme_bw() +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10",
        limits = c(1.0001, NA)
    ) +
    ylab("Run Accession") +
    ggtitle("Box Plot of Run Accession-NTpI after Normalization")
ggsave(
    "scaled_run_level_abundance_lim.pdf",
    g,
    width = 8,
    height = 5
)

g <- ggplot(model_sorted_all_gep_data_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance,
        y = RunAccession,
        fill = SequencerModel
    )) +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10"
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all runs")
ggsave("gep_density_run.pdf", g, width = 8, height = 5)

g <- ggplot(model_sorted_all_gep_data_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance,
        y = RunAccession,
        fill = SequencerModel
    )) +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10",
        limits = c(1.0001, NA)
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all runs")
ggsave("gep_density_run_lim.pdf", g, width = 8, height = 5)


g <- ggplot(all_gep_data_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance,
        y = SequencerModel,
        fill = SequencerModel
    )) +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10"
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all sequencers")
ggsave("gep_density_sequencer.pdf", g, width = 8, height = 5)


g <- ggplot(all_gep_data_normalized_long) +
    geom_density_ridges_gradient(aes(
        x = Abundance,
        y = SequencerModel,
        fill = SequencerModel
    )) +
    scale_x_continuous(
        name = "Normalized NTpI, log10 transformed",
        trans = "log10",
        limits = c(1.0001, NA)
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Normalized depth density distribution accross all sequencers")
ggsave("gep_density_sequencer_lim.pdf", g, width = 8, height = 5)


#' Sample Distances and Clustering
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
                all_gep_data_normalized[[accession1]],
                all_gep_data_normalized[[accession2]]
            ),
            method = "euclidean"
        )
    }
}

pheatmap(
    sample_distance_matrix,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename = "run_level_distance.pdf",
    width = 10,
    height = 8,
    main = "Run-level GEP difference (Euclidean)"
)

all_gep_data_sequencer <- all_gep_data_normalized

for (model_name in all_sequencer) {
    accessions <- datasets_metadata %>%
        dplyr::filter(SequencerModel == model_name) %>%
        dplyr::select(RunAccession) %>%
        unlist()
    all_gep_data_sequencer <- all_gep_data_sequencer %>%
        dplyr::mutate(!!model_name := rowMeans(select(., all_of(accessions)))) %>%
        dplyr::select(!all_of(accessions))
    rm(accessions)
}

sequencer_distance_matrix <- matrix(nrow = length(all_sequencer), ncol = length(all_sequencer))
colnames(sequencer_distance_matrix) <- all_sequencer
rownames(sequencer_distance_matrix) <- all_sequencer
for (accession1 in all_sequencer) {
    for (accession2 in all_sequencer) {
        sequencer_distance_matrix[accession1, accession2] <- dist(
            rbind(
                all_gep_data_sequencer[[accession1]],
                all_gep_data_sequencer[[accession2]]
            ),
            method = "euclidean"
        )
    }
}

pdf(
    "sequencer_level_distance.pdf",
    width = 10,
    height = 8
)
corrplot(
    sequencer_distance_matrix,
    is.corr = FALSE,
    type = "upper",
    title = "Sequencer-level GEP difference (Euclidean)"
)
dev.off()

#' Fit Zipf's law
#' Make data headless, which is also used for aggregation
#' and sorting.
all_gep_data_sorted <- all_gep_data_normalized
for (i in 10:length(all_gep_data_normalized)) {
    all_gep_data_sorted[, i] <- sort(as.double(unlist(all_gep_data_sorted[, i])))
}
all_gep_data_zipf_depth_long <- all_gep_data_sorted %>%
    dplyr::select(!(TRANSCRIPT_ID:GC)) %>%
    dplyr::mutate(RANK = n_of_records:1) %>%
    dplyr::sample_n(5000) %>%
    tidyr::gather(
        key = "RunAccession",
        value = "Abundance",
        -RANK
    )

g <- ggplot(all_gep_data_zipf_depth_long) +
    geom_point(aes(
        y = Abundance,
        x = RANK,
        color = RunAccession
    ), size = 0.2) +
    geom_function(
        data = NULL,
        fun = function(x) { 100000 / x }
    ) +
    scale_x_continuous(
        name = "log(Rank)",
        trans = "log10",
        labels = label_number()
    ) +
    scale_y_continuous(
        name = "log(Normalized NIpT)",
        trans = "log10",
        labels = label_number()
    ) +
    theme_bw() +
    ggtitle("Normalized depth for Zipf's Law accross all samples (Downsampled to 5000 each run)")
ggsave("gep_zipf.pdf", g, width = 8, height = 10)

g <- ggplot(all_gep_data_zipf_depth_long) +
    geom_point(aes(
        y = Abundance,
        x = RANK
    ), size = 0.2) +
    geom_function(
        data = NULL,
        fun = function(x) { 100000 / x }
    ) +
    scale_x_continuous(
        name = "log(Rank)",
        trans = "log10",
        labels = label_number()
    ) +
    scale_y_continuous(
        name = "log(Normalized NIpT)",
        trans = "log10",
        labels = label_number()
    ) +
    theme_bw() +
    facet_wrap(RunAccession ~ .) +
    ggtitle("Normalized depth for Zipf's Law accross all samples (Downsampled to 5000 each run)")
ggsave("gep_zipf.pdf", g, width = 8, height = 10)
