load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")

yasim_tsv_col_types <- cols(
gene_name=col_character(),
depth=col_double()
)
salmon_quant_sf_col_types <- cols(
Name=col_character(),
Length=col_double(),
EffectiveLength=col_double(),
TPM=col_double(),
NumReads=col_double()
)


