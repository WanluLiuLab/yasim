library("tidyverse")
library("ggridges")
library("arrow")


if (file.exists("all_fastq_data_sampled.parquet")) {
    all_data <- arrow::read_parquet("all_fastq_data_sampled.parquet")
} else {
    fns <- Sys.glob("ce11_as_2*.fq.stats.d")
    conditions <- fns %>%
        stringr::str_replace("ce11_", "") %>%
        stringr::str_replace(".fq.stats.d", "")

    all_data <- NULL
    for (i in seq_along(fns)) {
        this_data <- readr::read_tsv(
            file.path(fns[i], "all.tsv"),
            col_types = c(
                SEQID = col_character(),
                GC = col_double(),
                LEN = col_integer(),
                MEANQUAL = col_double()
            ),
            progress = FALSE,
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
        rm(i, this_data)
        gc()
    }

    as_ref <- readr::read_tsv(
        "../ce11_as_2.gtf.transcripts.tsv"
    )

    all_data_joint <- all_data %>%
        dplyr::mutate(
            TRANSCRIPT_ID = sapply(strsplit(all_data$SEQID, ":"), function(s) { s[1] })
        ) %>%
        dplyr::inner_join(as_ref, by = "TRANSCRIPT_ID") %>%
        dplyr::mutate(
            READ_COMPLETENESS = LEN / TRANSCRIBED_LENGTH
        )
    arrow::write_parquet(all_data_joint, "all_fastq_data_sampled.parquet")
}


g <- ggplot(all_data_joint) +
    geom_density_ridges_gradient(
        aes(
            x = READ_COMPLETENESS,
            y = Condition
        )
    ) +
    xlim(c(0, 1.5)) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Read Complexity of all conditions")

ggsave("fastq_rc_all.pdf", g, width = 8, height = 10)
