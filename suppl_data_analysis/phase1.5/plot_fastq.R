library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob("ce11_*.fq.stats.d")
conditions <- fns %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq.stats.d", "")

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
        dplyr::mutate(
            Condition = conditions[i]
        )
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
}
transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.chr1.gtf.transcripts.tsv",
    show_col_types = FALSE
)

all_data_merged <- all_data %>%
    tidyr::separate(
        "SEQID",
        c("TRANSCRIPT_ID", "READ_ID", "INPUT_DEPTH", "TmpCondition"),
        sep = ":"
    ) %>%
    dplyr::select(!(TmpCondition)) %>%
    dplyr::rename(READ_GC = GC) %>%
    dplyr::mutate(READ_ID = as.integer(READ_ID)) %>%
    dplyr::mutate(INPUT_DEPTH = as.double(INPUT_DEPTH)) %>%
    dplyr::inner_join(transcript_stats, by = "TRANSCRIPT_ID") %>%
    dplyr::mutate(READ_COMPLETENESS = LEN / TRANSCRIBED_LENGTH)


conditions <- unique(all_data_merged$Condition)

arrow::write_parquet(all_data_merged, "all_fastq_data.parquet")

g <- ggplot(all_data_merged) +
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

g <- ggplot(all_data_merged) +
    geom_density_ridges_gradient(
        aes(
            x = READ_GC,
            y = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("GC of all conditions")

ggsave("fastq_gc_all.pdf", g, width = 8, height = 5)

g <- ggplot(all_data_merged) +
    geom_density_ridges_gradient(
        aes(
            x = MEANQUAL,
            y = Condition
        )
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Mean Read Quality of all conditions")

ggsave("fastq_qual_all.pdf", g, width = 8, height = 5)
