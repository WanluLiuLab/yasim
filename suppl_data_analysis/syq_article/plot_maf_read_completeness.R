library("tidyverse")
library("ggridges")
library("arrow")

all_data <- NULL

fns <- Sys.glob("real_stats/*.fastq_trans.maf.gz.rlen.tsv")
conditions <- fns %>%
    stringr::str_replace("real_stats/", "") %>%
    stringr::str_replace(".fastq_trans.maf.gz.rlen.tsv", "")

metadata <- readr::read_csv(
    "metadata.csv",
    col_types = c(
        SequencerManufacturer = col_character(),
        SequencerModel = col_character(),
        Species = col_character(),
        SampleName = col_character(),
        Paper = col_character(),
        Mode = col_character(),
        Chemistry = col_character(),
        Basecaller = col_character(),
        Depth = col_double()
    ),
    comment = "#"
)

all_transcript_stats <- NULL

all_transcript_stats_fns <- Sys.glob("gtf/*")
species <- all_transcript_stats_fns %>%
    stringr::str_replace("gtf/", "") %>%
    stringr::str_replace(".tsv", "")

for (i in seq_along(all_transcript_stats_fns)) {
    this_data <- readr::read_tsv(
        all_transcript_stats_fns[i],
        col_types = c(
            TRANSCRIPT_ID = col_character(),
            GENE_ID = col_character(),
            NAIVE_LENGTH = col_integer(),
            NAIVE_LENGTH = col_integer(),
            EXON_NUMBER = col_integer()
        ),
        progress = TRUE,
        quote = "\'"
    ) %>%
        dplyr::mutate(
            Species = species[i]
        )
    if (is.null(all_transcript_stats)) {
        all_transcript_stats <- this_data
    } else {
        all_transcript_stats <- all_transcript_stats %>%
            dplyr::rows_append(this_data)
    }
    message(sprintf("Processing %d/%d", i, length(all_transcript_stats_fns)))
    rm(this_data, i)
    gc()
}

for (i in seq_along(fns)) {
    message(sprintf("Processing %d/%d", i, length(conditions)))
    this_data <- readr::read_tsv(
        fns[i],
        col_types = c(
            ALIGNED_TRANSCRIPT_ID = col_character(),
            READ_LENGTH = col_integer()
        ),
        progress = TRUE,
        quote = "\'"
    ) %>%
        dplyr::mutate(
            Condition = conditions[i]
        ) %>%
        dplyr::inner_join(
            metadata,
            by = c("Condition" = "SampleName")
        ) %>%
        dplyr::inner_join(
            all_transcript_stats,
            by = c("Species" = "Species", "ALIGNED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
        )
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
    rm(this_data, i)
    gc()
}
arrow::write_parquet(all_data, "all_maf_rc_data.parquet")

all_data_mutated <- all_data %>%
    dplyr::transmute(
        Condition = Condition,
        TRANSCRIPT_ID = ALIGNED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )

g <- ggplot(all_data_mutated) +
    geom_density_ridges(
        aes(
            x = READ_COMPLETENESS,
            y = Condition
        ),
        alpha = 0
    ) +
    ylab("density") +
    xlim(c(0.7, 1)) +
    theme_ridges() +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness.pdf", g, width = 10, height = 8)

g <- ggplot(all_data_mutated) +
    geom_histogram(
        aes(
            x = READ_COMPLETENESS
        )
    ) +
    ylab("density") +
    facet_wrap(. ~ Condition, scales = "free") +
    theme_ridges() +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness_hist.pdf", g, width = 20, height = 16)

