library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob("ce11_*.fq.bam.stats.d")
conditions <- fns %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq.bam.stats.d", "")

for (i in seq_along(conditions)) {
    if (is.null(all_data)) {
        all_data <- readr::read_tsv(
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
    }
    all_data <- all_data %>%
        dplyr::rows_append(
            readr::read_tsv(
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
        )
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

ggsave("sam_map_stat_all.pdf", g, width = 12, height = 8)

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = MAPPING_QUALITY,
            y = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("MAPQ of all conditions")

ggsave("sam_mapq_all.pdf", g, width = 12, height = 8)
