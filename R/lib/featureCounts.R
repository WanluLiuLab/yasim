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

#' Get featureCounts data and parse them into standard form
get_featureCounts_data <- function(featureCounts_tsv, n, ATTR_NAME) {
    message(sprintf("Reading %s...", featureCounts_tsv))
    featureCounts_data <- read_tsv(
        featureCounts_tsv,
        quote = "\'",
        col_types = featureCounts_tsv_col_types,
        col_names = featureCounts_tsv_col_names,
        comment = "#"
    )
    featureCounts_data <- featureCounts_data %>%
        dplyr::transmute(
            !!rlang::sym(ATTR_NAME) := Geneid,
            !!rlang::sym(sprintf("FEATURECOUNTS_%s_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", featureCounts_tsv))
    return(featureCounts_data)
}
