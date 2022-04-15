load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")


depth_data_col_type <- cols(
    TRANSCRIPT_ID=col_character(),
    BASE=col_number(),
    DEPTH=col_number()
)

depth_data_col_name <- c(
    "TRANSCRIPT_ID",
    "BASE",
    "DEPTH"
)


yasim_fa_stats_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    GENE_ID = col_character(),
    SEQNAME = col_character(),
    START = col_number(),
    END = col_number(),
    STRAND = col_character(),
    LEN = col_number(),
    GC = col_number()
)
