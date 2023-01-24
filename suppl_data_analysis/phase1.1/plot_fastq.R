library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- c(
    "ce11_badread_nanopore2018.fq.stats.d",
    "ce11_badread_nanopore2020.fq.stats.d",
    "ce11_badread_pacbio2016.fq.stats.d",
    "ce11_pbsim2_r94.fq.stats.d",
    "ce11_pbsim_clr.fq.stats.d"
)
conditions <- c(
    "NANOPORE2018",
    "NANOPORE2020",
    "PACBIO2016",
    "R94",
    "CLR"
)

for (i in seq_along(conditions)) {
    if (is.null(all_data)) {
        all_data <- readr::read_tsv(
            file.path(fns[i], "all.tsv"),
            col_types = c(
                SEQID = col_character(),
                GC = col_double(),
                LEN = col_integer(),
                MEANQUAL = col_double()
            ),
            progress = TRUE
        ) %>%
            dplyr::select(!(SEQID)) %>%
            dplyr::mutate(
                Condition = conditions[i]
            )
    }
    all_data <- all_data %>%
        dplyr::rows_append(
            readr::read_tsv(
                file.path(fns[i], "all.tsv"),
                col_types = c(
                    SEQID = col_character(),
                    GC = col_double(),
                    LEN = col_integer(),
                    MEANQUAL = col_double()
                ),
                progress = TRUE
            ) %>%
                dplyr::select(!(SEQID)) %>%
                dplyr::mutate(
                    Condition = conditions[i]
                )
        )
}

arrow::write_parquet(all_data, "all_fastq_data.parquet")

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = LEN,
            y = Condition,
            fill = Condition
        )
    ) +
    scale_x_continuous(
        name = "Length, log10 transformed",
        trans = "log10",
        limits = c(10, 10000)
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Length of all conditions")

ggsave("fastq_length_all.pdf", g, width = 12, height = 8)

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = GC,
            y = Condition,
            fill = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("GC of all conditions")

ggsave("fastq_gc_all.pdf", g, width = 12, height = 8)

g <- ggplot(all_data) +
    geom_density_ridges_gradient(
        aes(
            x = MEANQUAL,
            y = Condition,
            fill = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Mean Read Quality of all conditions")

ggsave("fastq_qual_all.pdf", g, width = 12, height = 8)
