library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob(file.path("real_stats", "*.fastq.sam.stats.d"))
conditions <- fns %>%
    stringr::str_replace("real_stats/", "") %>%
    stringr::str_replace(".fastq.sam.stats.d", "")

for (i in seq_along(conditions)) {
    this_data <- readr::read_tsv(
        file.path(fns[i], "read_stat.tsv"),
        col_types = c(
            QUERY_NAME = col_character(),
            MAP_STAT = col_character(),
            QUERY_LENGTH = col_integer(),
            REFERENCE_LENGTH = col_integer(),
            CIGAR_INFERRED_QUERY_LENGTH = col_integer(),
            CIGAR_INFERRED_READ_LENGTH = col_integer(),
            MAPPING_QUALITY = col_double()
        ),
        progress = TRUE
    ) %>%
        dplyr::select(!(QUERY_NAME)) %>%
        dplyr::mutate(
            Condition = conditions[i]
        )
    if (is.null(all_data)) {
        all_data <- this_data
    }
    all_data <- all_data %>%
        dplyr::rows_append(this_data)
    rm(this_data, i)
}

arrow::write_parquet(all_data, "all_sam_data.parquet")

g <- ggplot(all_data) +
    geom_bar(
        aes(
            y = Condition,
            fill = MAP_STAT
        ),
        stat = "count"
    ) +
    ylab("density") +
    theme_bw() +
    ggtitle("Mapping Status of all conditions")

ggsave("sam_map_stat_all.pdf", g, width = 8, height = 5)

g <- ggplot(all_data) +
    geom_bar(
        aes(
            y = Condition,
            fill = MAP_STAT
        ),
        stat = "count",
        position = "fill"
    ) +
    ylab("density") +
    theme_bw() +
    ggtitle("Mapping Status of all conditions")

ggsave("sam_map_stat_all_fill.pdf", g, width = 8, height = 5)
