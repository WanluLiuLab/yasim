cpptetgs_tsv_col_types <- cols(
    Name = col_character(),
    NumReads = col_double()
)
cpptetgs_tsv_col_names <- c(
    "Name",
    "NumReads"
)

get_cpptetgs_data <- function(cpptetgs_tsv, n, ATTR_NAME) {
    message(sprintf("Reading %s...", cpptetgs_tsv))
    cpptetgs_data <- read_tsv(
        cpptetgs_tsv,
        col_types = cpptetgs_tsv_col_types,
        col_names = cpptetgs_tsv_col_names,
        comment = "#"
    ) %>%
        dplyr::transmute(
            !!rlang::sym(ATTR_NAME) := Name,
            !!rlang::sym(sprintf("CPPTETGS_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", cpptetgs_tsv))
    return(cpptetgs_data)
}
