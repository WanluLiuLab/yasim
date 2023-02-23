library("tidyverse")
library("ggridges")
library("arrow")


if (file.exists("all_fastq_data_sampled.parquet")){
    all_data <- arrow::read_parquet("all_fastq_data_sampled.parquet")
} else{
    fns <- Sys.glob("ce11_*.fq.gz.stats.d")
    conditions <- fns %>%
        stringr::str_replace("ce11_", "") %>%
        stringr::str_replace(".fq.gz.stats.d", "")

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
            ) %>%
            dplyr::select(!(SEQID)) %>%
            dplyr::sample_n(10000)
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

    arrow::write_parquet(all_data, "all_fastq_data_sampled.parquet")
}


g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = LEN,
            y = Condition
        )
    ) +
    xlim(c(0, 3000)) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Length of all conditions")

ggsave("fastq_length_all.pdf", g, width = 8, height = 10)

g <- ggplot(all_data) +
    geom_boxplot(
        aes(
            x = LEN,
            y = Condition
        )
    ) +
    xlim(c(0, 3000)) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Length of all conditions")

ggsave("fastq_length_all_box.pdf", g, width = 8, height = 10)

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = GC,
            y = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("GC of all conditions")

ggsave("fastq_gc_all.pdf", g, width = 8, height = 5)

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = MEANQUAL,
            y = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Mean Read Quality of all conditions")

ggsave("fastq_qual_all.pdf", g, width = 8, height = 5)
