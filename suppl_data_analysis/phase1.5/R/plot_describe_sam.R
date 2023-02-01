library(tidyverse)
argv <- commandArgs(trailingOnly = TRUE)

outdir_path <- argv[1]

read_stat_tsv <- readr::read_tsv(
    file.path(outdir_path, "read_stat.tsv"),
    col_types = cols(
        QUERY_NAME = col_character(),
        MAP_STAT = col_character(),
        QUERY_LENGTH = col_integer(),
        REFERENCE_LENGTH = col_integer(),
        CIGAR_INFERRED_QUERY_LENGTH = col_integer(),
        CIGAR_INFERRED_READ_LENGTH = col_integer(),
        MAPPING_QUALITY = col_integer(),
    ),
    progress = TRUE
)

map_stat_plot <- ggplot(read_stat_tsv) +
    geom_bar(aes(y = MAP_STAT), stat = "count") +
    theme_bw() +
    ggtitle("Alignment plot")
ggsave(
    file.path(outdir_path, "map_stat_plot.pdf"),
    map_stat_plot,
    height = 8,
    width = 10
)

primiary_read_stat_tsv <- read_stat_tsv %>%
    dplyr::filter(
        MAP_STAT == "primiary"
    )

length_distribution <- primiary_read_stat_tsv %>%
    tidyr::gather(
        key = "length_type",
        value = "value",
        -QUERY_NAME,
        -MAP_STAT,
        -MAPPING_QUALITY
    )

length_distributions_plot <- ggplot(length_distribution) +
    geom_density(aes(x = value, color = length_type)) +
    theme_bw() +
    xlim(c(0, 5000)) +
    ggtitle("Length distribution plot in primiary alignments (x limited to 5000)")
ggsave(
    file.path(outdir_path, "length_distributions_plot.pdf"),
    length_distributions_plot,
    height = 8,
    width = 10
)

pileup_stat_tsv <- readr::read_tsv(
    file.path(outdir_path, "pileup_stat.tsv.gz"),
    col_types = cols(
        REFERENCE_NAME = col_character(),
        REFERENCE_POS = col_double(),
        NUM_READS = col_double()
    ),
    progress = TRUE,
    n_max = 10000000
) %>%
    dplyr::filter(NUM_READS > 0)

mean_depth <- mean(pileup_stat_tsv$NUM_READS)

depth_plot <- ggplot(pileup_stat_tsv) +
    geom_histogram(aes(x = NUM_READS), binwidth = 1) +
    theme_bw() +
    xlim(c(0, 2 * mean_depth)) +
    ggtitle(
        sprintf(
            "per-base Depth distribution (data limited to 10000000, x limited to 2*mean), mean %f",
            mean(pileup_stat_tsv$NUM_READS)
        )
    )
ggsave(
    file.path(outdir_path, "depth_plot.pdf"),
    depth_plot,
    height = 8,
    width = 10
)
