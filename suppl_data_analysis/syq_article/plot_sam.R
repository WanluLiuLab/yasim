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

g <- metadata %>%
    dplyr::mutate(Paper = factor(
        Paper,
        level = sort(unique(.$Paper))
    )) %>%
    dplyr::arrange(desc(Paper)) %>%
    dplyr::mutate(SampleName = factor(
        SampleName,
        level = unique(.$SampleName)
    )) %>%
    ggplot() +
    geom_bar(
        aes(
            x = Depth,
            y = SampleName,
            fill = Paper
        ),
        stat = "identity"
    ) +
    theme_ridges() +
    ggtitle("Sequencing Depth of all conditions")

ggsave("sam_depth_all.pdf", g, width = 10, height = 8)

for (i in seq_along(conditions)) {
    message(sprintf("Reading %s --  %d/%d", fns[i], i, length(conditions)))
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
        na = "None" # Be compatible with lower versions of labw_utils
    ) %>%
        dplyr::select(!(QUERY_NAME)) %>%
        dplyr::mutate(across(where(is.numeric), replace_na, 0)) %>%
        dplyr::mutate(
            Condition = conditions[i]
        )
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
    rm(this_data, i)
}
all_data <- all_data %>%
    dplyr::inner_join(metadata, by = c("Condition" = "SampleName"))

arrow::write_parquet(all_data, "all_sam_data.parquet")

all_data_sampled <- all_data %>%
    dplyr::sample_frac(0.01)

g <- ggplot(all_data_sampled) +
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

ggsave("sam_map_stat_all.pdf", g, width = 8, height = 10)

g <- ggplot(all_data_sampled) +
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

ggsave("sam_map_stat_all_fill.pdf", g, width = 8, height = 10)

all_data_sampled_primiary <- all_data_sampled %>%
    dplyr::filter(MAP_STAT == "primiary")

g <- ggplot(all_data_sampled_primiary) +
    geom_hex(aes(x = QUERY_LENGTH, y = REFERENCE_LENGTH)) +
    theme_bw() +
    scale_fill_continuous(trans = "log10") +
    facet_wrap(. ~ Condition)
ggsave("sam_query_reference_relation.pdf", g, width = 10, height = 8)

