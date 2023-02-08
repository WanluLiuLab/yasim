library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob("ce11_*.fq.gz_trans.maf.gz.rlen.tsv")
conditions <- fns %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq.gz_trans.maf.gz.rlen.tsv", "")

for (i in seq_along(fns)) {
    this_data <- readr::read_tsv(
        fns[i],
        col_types = c(
            ALIGNED_TRANSCRIPT_ID = col_character(),
            SIMULATED_TRANSCRIPT_ID = col_character(),
            READ_LENGTH = col_integer()
        ),
        progress = TRUE,
        quote = "\'"
    ) %>%
        dplyr::mutate(
            Condition = conditions[i]
        )
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
    message(sprintf("Processing %d/%d", i, length(conditions)))
    rm(this_data, i)
    gc()
}
aligned_transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.chr1.gtf.transcripts.tsv",
    show_col_types = FALSE
)
simulated_transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq_as.chr1.gtf.transcripts.tsv",
    show_col_types = FALSE
)
arrow::write_parquet(all_data, "all_last_aligned_maf_data.parquet")

correctness <- all_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
        correctness = sum(SIMULATED_TRANSCRIPT_ID == ALIGNED_TRANSCRIPT_ID) / n()
    ) %>%
    dplyr::ungroup()

correctness_mutated <- correctness %>%
    tidyr::separate(
        "Condition",
        c("SIMULATOR", "MODE", "DGEID", "DIUID", "REPID")
    )

g <- ggplot(correctness_mutated) +
    geom_boxplot(
        aes(
            y = sprintf("%s_%s", SIMULATOR, MODE),
            x = correctness
        )
    ) +
    xlim(c(0.0, 1.0)) +
    theme_bw() +
    facet_wrap(DGEID ~ DIUID) +
    ggtitle("Primiary Mapping Rate of all conditions")
ggsave("last_correctness.pdf", g, width = 5, height = 4)

all_data_sim <- all_data %>%
    dplyr::inner_join(
        simulated_transcript_stats,
        by = c("SIMULATED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::sample_n(10000) %>%
    dplyr::transmute(
        Condition = sprintf("%s_sim", Condition),
        TRANSCRIPT_ID = SIMULATED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )
all_data_aln <- all_data %>%
    dplyr::inner_join(
        aligned_transcript_stats,
        by = c("ALIGNED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::sample_n(10000) %>%
    dplyr::transmute(
        Condition = sprintf("%s_aln", Condition),
        TRANSCRIPT_ID = ALIGNED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )

all_data_long <- dplyr::rows_append(all_data_sim, all_data_aln)
rm(all_data_sim, all_data_aln)
gc()

all_data_long_mutated <- all_data_long %>%
    tidyr::separate(
        "Condition",
        c("SIMULATOR", "MODE", "DGEID", "DIUID", "REPID", "SOURCE")
    )

g <- ggplot(all_data_long_mutated) +
    geom_density_ridges(
        aes(
            x = READ_COMPLETENESS,
            y = sprintf("%s_%s", SIMULATOR, MODE),
            color = SOURCE
        ),
        alpha = 0
    ) +
    ylab("density") +
    xlim(c(0.7, 1)) +
    theme_ridges() +
    facet_wrap(DGEID ~ DIUID) +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness.pdf", g, width = 8, height = 8)
