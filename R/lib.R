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
cpptetgs_tsv_col_types <- cols(
    Name = col_character(),
    NumReads = col_double()
)
cpptetgs_tsv_col_names <- c(
    "Name",
    "NumReads"
)
stringtie_quant_tsv_col_types <- cols(
    gene_id = col_character(),
    transcript_id = col_character(),
    reference_id = col_character(),
    ref_gene_id = col_character(),
    ref_transcript_id = col_character(),
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
yasim_unmapped_stats_col_types <- cols(
    TRANSCRIPT_ID=col_character(),
    SUCCESS_ALN=col_number(),
    FAILED_ALN=col_number()
)
htseq_count_tsv_col_types <- cols(
    Name = col_character(),
    NumReads = col_double()
)
htseq_count_tsv_col_names <- c(
    "Name",
    "NumReads"
)
ss_tsv_col_types <- cols(
    LEN = col_number(),
    GC = col_number(),
    IS_MAPPED = col_number(),
    MAPQ = col_number(),
    ALNQ = col_number()
)

#' Get featureCounts data and parse them into standard form
get_featureCounts_data <- function(featureCounts_tsv, n) {
    message(sprintf("Reading %s...", featureCounts_tsv))
    featureCounts_data <- read_tsv(
        featureCounts_tsv,
        quote = "\'",
        col_types = featureCounts_tsv_col_types,
        col_names = featureCounts_tsv_col_names,
        comment = "#"
    )
    featureCounts_data <- featureCounts_data %>%
        dplyr::filter(NumReads > 0) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = Geneid,
            !!rlang::sym(sprintf("FEATURECOUNTS_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", featureCounts_tsv))
    return(featureCounts_data)
}

get_fq_stats_data <- function(fq_stats, n) {
    message(sprintf("Reading %s...", fq_stats))
    fq_stats_data <- read_tsv(
        fq_stats,
        col_types = yasim_fq_stats_col_types,
        comment = "#"
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("YASIM_%d_INPUT_DEPTH", n)) := INPUT_DEPTH,
            !!rlang::sym(sprintf("YASIM_%d_ACTUAL_N_OF_READS", n)) := SIMULATED_N_OF_READS
        )
    message(sprintf("Reading %s... DONE", fq_stats))
    return(fq_stats_data)
}
get_unmapped_stats_data <- function(unmapped_stats, n) {
    message(sprintf("Reading %s...", unmapped_stats))
    unmapped_stats_data <- read_tsv(
        unmapped_stats,
        col_types = yasim_unmapped_stats_col_types,
        comment = "#"
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("YASIM_UNMAPPED_%d_ACTUAL_N_OF_READS", n)) := FAILED_ALN,
            !!rlang::sym(sprintf("YASIM_MAPPED_%d_ACTUAL_N_OF_READS", n)) := SUCCESS_ALN
        )
    message(sprintf("Reading %s... DONE", fq_stats))
    return(unmapped_stats_data)
}

get_salmon_data <- function(salmon_quant_sf, n) {
    message(sprintf("Reading %s...", salmon_quant_sf))
    salmon_quant_sf_data <- read_tsv(salmon_quant_sf, col_types = salmon_quant_sf_col_types) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = Name,
            !!rlang::sym(sprintf("SALMON_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", salmon_quant_sf))
    return(salmon_quant_sf_data)
}

get_cpptetgs_data <- function(cpptetgs_tsv, n) {
    message(sprintf("Reading %s...", cpptetgs_tsv))
    cpptetgs_data <- read_tsv(
        cpptetgs_tsv,
        col_types = cpptetgs_tsv_col_types,
        col_names = cpptetgs_tsv_col_names,
        comment = "#"
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = Name,
            !!rlang::sym(sprintf("CPPTETGS_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", cpptetgs_tsv))
    return(cpptetgs_data)
}



get_htseq_count_data <- function(htseq_count_tsv, n) {
    htseq_count_data <- read_tsv(
        htseq_count_tsv,
        quote = "\'",
        col_types = htseq_count_tsv_col_types,
        col_names = htseq_count_tsv_col_names
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = Name,
            !!rlang::sym(sprintf("HTSEQ_COUNT_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", htseq_count_tsv))
    return(htseq_count_data)
}
