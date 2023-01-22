library(tidyverse)


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

yasim_fa_stats_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    GENE_ID = col_character(),
    SEQNAME = col_character(),
    START = col_number(),
    END = col_number(),
    STRAND = col_character(),
    ABSOLUTE_LENGTH = col_number(),
    TRANSCRIBED_LENGTH = col_number(),
    GC = col_number()
)

nano_depth_data <- read_tsv(
    "NANO_TRANS.bam.depth.tsv",
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(NANOPORE_AVG_DEPTH = mean(DEPTH))
pacb_depth_data <- read_tsv(
    "PACB_TRANS.bam.depth.tsv",
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(PACB_AVG_DEPTH = mean(DEPTH))
illm_depth_data <- read_tsv(
    "ILLM_TRANS.bam.depth.tsv",
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
) %>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(ILLM_AVG_DEPTH = mean(DEPTH))


fa_stats_data <- read_tsv(
    "ce11.reference_transcripts.fa.stats",
    col_types = yasim_fa_stats_col_types
)

all_data <- fa_stats_data %>%
    dplyr::full_join(pacb_depth_data, by = "TRANSCRIPT_ID") %>%
    dplyr::full_join(nano_depth_data, by = "TRANSCRIPT_ID") %>%
    dplyr::full_join(illm_depth_data, by = "TRANSCRIPT_ID") %>%
    dplyr::mutate(across(where(is.numeric), replace_na, 0))

write_tsv(all_data, "all_data.tsv")
