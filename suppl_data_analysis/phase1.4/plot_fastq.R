library("tidyverse")
library("pheatmap")
library("scales")
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
    rm(i,this_data)
    gc()
}

arrow::write_parquet(all_data, "all_fastq_data.parquet")

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
