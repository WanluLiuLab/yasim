htseq_count_tsv_col_types <- cols(
    Name = col_character(),
    NumReads = col_double()
)
htseq_count_tsv_col_names <- c(
    "Name",
    "NumReads"
)

get_htseq_count_data <- function(htseq_count_tsv, n, ATTR_NAME) {
    htseq_count_data <- read_tsv(
        htseq_count_tsv,
        quote = "\'",
        col_types = htseq_count_tsv_col_types,
        col_names = htseq_count_tsv_col_names
    ) %>%
        dplyr::transmute(
            !!rlang::sym(ATTR_NAME) := Name,
            !!rlang::sym(sprintf("HTSEQ_COUNT_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", htseq_count_tsv))
    return(htseq_count_data)
}

