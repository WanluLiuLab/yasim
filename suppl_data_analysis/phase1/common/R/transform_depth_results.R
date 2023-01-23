file_description <- "transform_depth_results.R -- Transform per-base depth data into per-transcript."

#' This file is used to merge TSVs produced by Transcript-Level Quantifiers for Downstream DGE Analysis.
#' The

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--input", help = "Input depth tsv produced by Samtools", type = "character")
p <- add_argument(p, "--output", help = "Output basename", type = "character")
argv <- parse_args(p)

source(argv$libfile)

depth_data <- read_tsv(
    argv$input,
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(AVG_DEPTH = mean(DEPTH))

write_tsv(all_data, argv$output)
