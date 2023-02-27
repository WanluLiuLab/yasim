library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

fns_genome <- Sys.glob("ce11_*_*_dge?_diu?_?.fq.gz.bam.stats.d")
conditions_genome <- fns_genome %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq.gz.bam.stats.d", "") %>%
    sprintf("%s_genome", .)

fns_transcriptome <- Sys.glob("ce11_*_*_dge?_diu?_?_trans.fq.gz.bam.stats.d")
conditions_transcriptome <- fns_transcriptome %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace("_trans", "_transcriptome") %>%
    stringr::str_replace(".fq.gz.bam.stats.d", "")

conditions <- c(conditions_genome, conditions_transcriptome)
fns <- c(fns_genome, fns_transcriptome)
rm(conditions_genome, conditions_transcriptome, fns_genome, fns_transcriptome)

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

all_data_mutated <- all_data %>%
    tidyr::separate(
        "Condition",
        c("SIMULATOR", "SEQUENCER", "MODE", "DGEID", "DIUID", "REPID", "ALNTO")
    )

all_data_alignment_rate <- all_data_mutated %>%
    dplyr::group_by(SIMULATOR, MODE, SEQUENCER, DGEID, DIUID, REPID, ALNTO) %>%
    dplyr::summarise(PRIMIARY_ALN_RATE = sum(MAP_STAT == "primiary") / n()) %>%
    dplyr::ungroup()

g <- ggplot(all_data_alignment_rate) +
    geom_boxplot(
        aes(
            y = sprintf("%s_%s_%s", SIMULATOR, SEQUENCER, MODE),
            x = PRIMIARY_ALN_RATE,
            color = ALNTO
        )
    ) +
    theme_bw() +
    facet_wrap(DGEID ~ DIUID) +
    ggtitle("Primiary Mapping Rate of all conditions")

ggsave("sam_primiary_mapping_rate.pdf", g, width = 10, height = 8)
