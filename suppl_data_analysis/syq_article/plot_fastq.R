library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob(file.path("real_stats", "*.fastq.stats.d"))
conditions <- fns %>%
    stringr::str_replace("real_stats/", "") %>%
    stringr::str_replace(".fastq.stats.d", "")
for (i in seq_along(fns)) {
    message(sprintf("Reading %s --  %d/%d", fns[i], i, length(conditions)))
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

all_data <- all_data %>%
    dplyr::inner_join(metadata, by = c("Condition" = "SampleName"))

arrow::write_parquet(all_data, "all_fastq_data.parquet")


all_data_sampled <- all_data %>%
    dplyr::sample_frac(0.01)

g <- ggplot(all_data_sampled) +
    geom_density_ridges_gradient(
        aes(
            x = LEN,
            y = Condition,
            fill = SequencerModel
        )
    ) +
    xlim(c(0, 3000)) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Length of all conditions")

ggsave("fastq_length_all.pdf", g, width = 8, height = 10)

g <- ggplot(all_data_sampled) +
    geom_density_ridges_gradient(
        aes(
            x = GC,
            y = Condition,
            fill = Species
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("GC of all conditions")

ggsave("fastq_gc_all.pdf", g, width = 8, height = 10)

g <- ggplot(all_data_sampled) +
    geom_density_ridges_gradient(
        aes(
            x = MEANQUAL,
            y = Condition,
            fill = Paper
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Mean Read Quality of all conditions")

ggsave("fastq_qual_all.pdf", g, width = 8, height = 10)

g <- ggplot(all_data_sampled) +
    geom_hex(aes(x = LEN, y = MEANQUAL)) +
    theme_bw() +
    facet_wrap(. ~ Condition)
ggsave("fastq_len_qual_relation.pdf", g, width = 10, height = 8)

g <- ggplot(all_data_sampled) +
    geom_hex(aes(x = LEN, y = GC)) +
    theme_bw() +
    facet_wrap(. ~ Condition)
ggsave("fastq_len_gc_relation.pdf", g, width = 10, height = 8)

g <- ggplot(all_data_sampled) +
    geom_hex(aes(x = GC, y = MEANQUAL)) +
    theme_bw() +
    facet_wrap(. ~ Condition)
ggsave("fastq_gc_qual_relation.pdf", g, width = 10, height = 8)

