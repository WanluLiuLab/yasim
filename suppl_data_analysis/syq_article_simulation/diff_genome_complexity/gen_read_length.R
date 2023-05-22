library("tidyverse")
library("arrow")


if (file.exists("all_fastq_data_rlen.parquet")) {
    quit()
} else {
    fns <- Sys.glob("ce11_as_*.fq.stats.d")
    conditions <- fns %>%
        stringr::str_replace("ce11_", "") %>%
        stringr::str_replace(".fq.stats.d", "")

    all_data <- NULL
    for (i in seq_along(fns)) {
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
            dplyr::mutate(
                Condition = conditions[i]
            ) %>%
            dplyr::select(
                Condition,
                LEN
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
    arrow::write_parquet(all_data, "all_fastq_data_rlen.parquet")
}
