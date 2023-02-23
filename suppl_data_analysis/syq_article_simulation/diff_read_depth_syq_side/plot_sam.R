library("tidyverse")
library("arrow")

fns_genome <- Sys.glob("ce11_*.fq.gz.sam.stats.d")
conditions_genome <- fns_genome %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq.gz.sam.stats.d", "")


conditions <- c(conditions_genome)
fns <- c(fns_genome)
rm(conditions_genome, fns_genome)

all_data <- NULL
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
        progress = FALSE,
        lazy = TRUE
    ) %>%
        dplyr::select(MAP_STAT) %>%
        dplyr::mutate(
            Condition = conditions[i]
        ) %>%
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

arrow::write_parquet(all_data, "all_sam_data_sampled.parquet")

all_data_alignment_rate <- all_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(PRIMIARY_ALN_RATE = sum(MAP_STAT == "primiary") / n()) %>%
    dplyr::ungroup()

g <- ggplot(all_data_alignment_rate) +
    geom_boxplot(
        aes(
            y = Condition,
            x = PRIMIARY_ALN_RATE
        )
    ) +
    theme_bw() +
    facet_wrap(.~Condition) +
    ggtitle("Primiary Mapping Rate of all conditions")

ggsave("sam_primiary_mapping_rate.pdf", g, width = 10, height = 8)

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
