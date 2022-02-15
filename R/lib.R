load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")

yasim_tsv_col_types <- cols(
    gene_name = col_character(),
    depth = col_double()
)
salmon_quant_sf_col_types <- cols(
    Name = col_character(),
    Length = col_double(),
    EffectiveLength = col_double(),
    TPM = col_double(),
    NumReads = col_double()
)
stringtie_quant_tsv_col_types <- cols(
    gene_id = col_character(),
    transcript_id = col_character(),
    reference_id = col_character(),
    ref_gene_id = col_character(),
    ref_gene_name = col_character(),
    cov = col_double(),
    FPKM = col_double(),
    TPM = col_double()
)
featureCounts_tsv_col_types <- cols(
    Geneid = col_character(),
    Chr = col_character(),
    Start = col_character(),
    End = col_character(),
    Strand = col_character(),
    Length = col_number(),
    NumReads = col_double(),
)
featureCounts_tsv_col_names <- c(
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length",
    "NumReads"
)
htseq_quant_tsv_col_types <- cols(
    Name=col_character(),
    NumReads=col_double()
)
htseq_quant_tsv_col_names <- c(
    "Name",
    "NumReads"
)
