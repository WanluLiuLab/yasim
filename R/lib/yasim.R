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

yasim_unmapped_stats_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    SUCCESS_ALN = col_number(),
    FAILED_ALN = col_number()
)


yasim_pileup_stats_col_types <- cols(
    REFERENCE_NAME = col_character(),
    REFERENCE_POS = col_number(),
    NUM_READS = col_number()
)
yasim_pileup_stats_merged_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    AVG_DEPTH = col_number()
)

get_fq_stats_data <- function(fq_stats, n) {
    message(sprintf("Reading %s...", fq_stats))
    fq_stats_data <- read_tsv(
        fq_stats,
        col_types = yasim_fq_stats_col_types,
        comment = "#"
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("YASIM_%s_INPUT_DEPTH", n)) := INPUT_DEPTH,
            !!rlang::sym(sprintf("YASIM_%s_ACTUAL_N_OF_READS", n)) := SIMULATED_N_OF_READS
        )
    message(sprintf("Reading %s... DONE", fq_stats))
    return(fq_stats_data)
}

get_pileup_stats_merged_data <- function(pileup_stats_merged, n) {
    message(sprintf("Reading %s...", pileup_stats_merged))
    pileup_stats_merged_data <- read_tsv(
        pileup_stats_merged,
        col_types = yasim_pileup_stats_merged_col_types,
        comment = "#"
    ) %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(sprintf("PILEUP_STATS_YASIM_%s_AVG_DEPTH", n)) := AVG_DEPTH
        )
    message(sprintf("Reading %s... DONE", pileup_stats_merged))
    return(pileup_stats_merged_data)
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
