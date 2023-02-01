library(tidyverse)

argv <- commandArgs(trailingOnly = TRUE)

input_depth_tsv <- argv[1]
output_depth_tsv <- argv[2]

depth_data <- readr::read_tsv(
    input_depth_tsv
) %>%
    dplyr::rename(TRANSCRIPT_ID = REFERENCE_NAME) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(AVG_DEPTH = mean(NUM_READS))

readr::write_tsv(depth_data, output_depth_tsv)
