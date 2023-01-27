library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

fns <- Sys.glob(file.path("real_stats", "*.fastq.sam.stats.d"))
conditions <- fns %>%
    stringr::str_replace("real_stats/", "") %>%
    stringr::str_replace(".fastq.sam.stats.d", "")

all_data <- NULL
for (i in seq_along(conditions)) {
    for (fn in Sys.glob(file.path(fns[i], "pileup_stat.tsv.gz.sampled.parquet.d", "*.parquet"))) {
        this_data <- arrow::read_parquet(fn) %>%
            dplyr::mutate(
                Condition = conditions[i]
            ) %>%
            dplyr::select(NUM_READS, Condition)
        if (is.null(all_data)) {
            all_data <- this_data
        }
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
}

all_data_filetred <- all_data %>%
    dplyr::filter(NUM_READS != 0)

arrow::write_parquet(all_data, "all_sam_data_depth.parquet")

g <- ggplot(all_data) +
    geom_histogram(
        aes(
            x = NUM_READS + 1
        )
    ) +
    scale_x_continuous(
        limits = c(NA, 100),
        trans = "log10"
    ) +
    theme_ridges() +
    facet_wrap(Condition ~ ., scales = "free") +
    ggtitle("Depth distribution of all conditions")

ggsave("sam_depth.pdf", g, width = 12, height = 8)

all_data_filetred %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(MEAN_DEPTH = mean(NUM_READS))
