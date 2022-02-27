load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")

yasim_depth_tsv_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    INPUT_DEPTH = col_double()
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
yasim_fq_stats_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    INPUT_DEPTH = col_number(),
    SIMULATED_N_OF_READS = col_number()
)
yasim_ground_truth_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    GENE_ID = col_character(),
    SEQNAME = col_character(),
    START = col_number(),
    END = col_number(),
    STRAND = col_character(),
    LEN = col_number(),
    GC = col_number(),
    INPUT_DEPTH = col_number(),
    SIMULATED_N_OF_READS = col_number(),
    SIMULATED_RPM = col_number(),
    SIMULATED_RPK = col_number(),
    SIMULATED_RPKM = col_number(),
    SIMULATED_TPM = col_number()
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
    ref_TRANSCRIPT_ID = col_character(),
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
# htseq_quant_tsv_col_types <- cols(
#     Name = col_character(),
#     NumReads = col_double()
# )
# htseq_quant_tsv_col_names <- c(
#     "Name",
#     "NumReads"
# )
ss_tsv_col_types <- cols(
    LEN=col_number(),
    GC=col_number(),
    IS_MAPPED=col_number(),
    MAPQ=col_number(),
    ALNQ=col_number()
)
