library("tidyverse")
library("pheatmap")
library("ggridges")
library("arrow")

all_data <- NULL

fns <- Sys.glob("*.fastq.stats.d")
conditions <- fns %>%
    stringr::str_replace(".fastq.stats.d", "")
for (i in seq_along(fns)) {
    this_data <- readr::read_tsv(
        file.path(fns[i], "all.tsv"),
        col_types = c(
            SEQID = col_character(),
            GC = col_double(),
            LEN = col_integer(),
            MEANQUAL = col_double()
        ),
        progress = TRUE,
        quote = "\'"
    ) %>%
        dplyr::select(!(SEQID)) %>%
        dplyr::mutate(Condition = conditions[i])
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
    rm(i, this_data)
    gc()
}

arrow::write_parquet(all_data, "all_fastq_data_length.parquet")

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

ggsave("fastq_length_all.pdf", g, width = 8, height = 5)

all_error <- readr::read_tsv(
    "all_last_mapq.tsv",
    col_types = c(
        FILENAME = col_character(),
        INSERTION = col_double(),
        DELETION = col_double(),
        MATCH = col_double(),
        SUBSTITUTION = col_double()
    )
) %>%
    dplyr::mutate(
        SampleName = FILENAME %>%
            stringr::str_replace(
                ".fastq_trans.maf",
                ""
            )
    ) %>%
    dplyr::select(!(FILENAME))
all_error_long <- all_error %>%
    tidyr::gather(
        key = "EventType",
        value = "EventCount",
        -SampleName
    )

g <- ggplot(all_error_long) +
    geom_bar(
        aes(
            y = SampleName,
            x = EventCount,
            fill = EventType
        ),
        stat = "identity"
    ) +
    theme_bw()
ggsave("maf_error_rate.pdf", g, width = 8, height = 5)


g <- ggplot(all_error_long) +
    geom_bar(
        aes(
            y = SampleName,
            x = EventCount,
            fill = EventType
        ),
        stat = "identity",
        position = "fill"
    ) +
    theme_bw()
ggsave("maf_error_rate_fill.pdf", g, width = 8, height = 5)

all_error_accuracy <- all_error %>%
    dplyr::mutate(Accuracy = MATCH / (MATCH + INSERTION + DELETION + SUBSTITUTION))

g <- ggplot(all_error_accuracy) +
    geom_bar(
        aes(
            y = SampleName,
            x = Accuracy
        ),
        stat = "identity"
    ) +
    theme_bw()
ggsave("maf_accuracy.pdf", g, width = 8, height = 5)

g <- ggplot(all_error_accuracy) +
    geom_point(
        aes(
            y = SampleName,
            x = Accuracy
        )
    ) +
    theme_bw()
ggsave("maf_accuracy_log.pdf", g, width = 8, height = 5)
