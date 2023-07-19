library("tidyverse")
library("arrow")

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
        dplyr::mutate(
            Condition = conditions[i]
        )
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

arrow::write_parquet(all_data, "all_sam_data.parquet")
