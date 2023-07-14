library(argparser)

p <- arg_parser("DESC")
p <- add_argument(p, "--input", help = "Input depth tsv produced by Samtools", type = "character")
p <- add_argument(p, "--output", help = "Transformed Isoform-Level depth tsv", type = "character")
p <- add_argument(p, "--trans_stats", help = "Input GTF trans.stats", type = "character")


argv <- parse_args(p)

load_package <- function(name) {
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")


depth_data_col_type <- cols(
    TRANSCRIPT_ID = col_character(),
    BASE = col_number(),
    DEPTH = col_number()
)

depth_data_col_name <- c(
    "TRANSCRIPT_ID",
    "BASE",
    "DEPTH"
)


trans_stats <- read_tsv(
    argv$trans_stats,
    show_col_types = FALSE
) %>%
    dplyr::select(TRANSCRIBED_LENGTH, TRANSCRIPT_ID, GENE_ID)

depth_data <- read_tsv(
    argv$input,
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
) %>%
    dplyr::filter(DEPTH!=0) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(SUM_BASE = sum(DEPTH)) %>%
    dplyr::filter(SUM_BASE!=0) %>%
    dplyr::select(TRANSCRIPT_ID, SUM_BASE) %>%
    dplyr::inner_join(trans_stats, by="TRANSCRIPT_ID") %>%
    dplyr::mutate(DEPTH = SUM_BASE / TRANSCRIBED_LENGTH)

arrow::write_parquet(depth_data, argv$output)
print(mean(depth_data$DEPTH))

