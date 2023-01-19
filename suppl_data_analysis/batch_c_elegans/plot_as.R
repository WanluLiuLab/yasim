library("tidyverse")
library("pheatmap")

as_event_types <- c("alt_3prime", "alt_5prime", "exon_skip", "intron_retention", "mult_exon_skip", "mutex_exons")

#' Get dataset metadata
datasets_metadata <-
    readr::read_csv(
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
        dplyr::filter(Generation != "TGS") %>%
        dplyr::mutate(as_filename = sprintf(
            "./datasets/%s/%s_SPLADDER/all_as.tsv.aggregated.tsv",
            SequencerModel,
            RunAccession
        )) %>%
        dplyr::mutate(
            alt_3prime = NULL,
            alt_5prime = NULL,
            exon_skip = NULL,
            intron_retention = NULL,
            mult_exon_skip = NULL,
            mutex_exons = NULL
        )

#' Raw sequence run accessions
all_run_accessions <- unlist(datasets_metadata$RunAccession)

#' Read data
for (i in seq_len(nrow(datasets_metadata))) {
    data_tuple <- datasets_metadata[i,]
    if (file.exists(data_tuple$as_filename)) {
        x <- readr::read_tsv(
            data_tuple$as_filename,
            show_col_types = FALSE,
            col_types = c(
                GENE_ID = col_character(),
                alt_3prime = col_integer(),
                alt_5prime = col_integer(),
                exon_skip = col_integer(),
                intron_retention = col_integer(),
                mult_exon_skip = col_integer(),
                mutex_exons = col_integer(),
                sum = col_integer()
            )
        )
        for (as_event_type in as_event_types) {
            datasets_metadata[[as_event_type]][i] <- sum(x[[as_event_type]])
        }
        rm(x)
    } else {
        warning(sprintf("%s not found!", data_tuple$filename))
    }
    rm(data_tuple)
}

datasets_metadata_long <- datasets_metadata %>%
    dplyr::select(!(Generation:ProjectAccession)) %>%
    dplyr::select(!(as_filename)) %>%
    tidyr::gather(key = "ASEventType", value = "Number", -RunAccession)

g <- ggplot(datasets_metadata_long) +
    geom_bar(aes(x = RunAccession, y = Number, fill = ASEventType), stat = "identity", position = "fill") +
    ggtitle("Ratio between alternative splicing events") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave(
    "ratio_as_events.pdf",
    g,
    width = 8, height = 5
)


datasets_matrix <- datasets_metadata %>%
    dplyr::select(!(Generation:as_filename)) %>%
    as.matrix()

row.names(datasets_matrix) <- datasets_metadata$RunAccession
pheatmap(
    log10(datasets_matrix + 1),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    filename = "count_as_events_heatmap.pdf",
    width = 10,
    height = 8,
    main = "Absolute count of AS Events (log10 +1 transformed)"
)
