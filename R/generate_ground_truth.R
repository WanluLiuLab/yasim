file_description <- "Generates ground truth from *.fq.stats and *.fa.stats"

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--fq_stats", help = "*.fq.stats produced by YASIM simulator", type = "character")
p <- add_argument(p, "--fa_stats", help = "*.fa.stats produced by YASIM transcriptor", type = "character")
p <- add_argument(p, "--output", help = "Output basename.", type = "character")
argv <- parse_args(p)

source(argv$libfile)

fq_stats_data <- read_tsv(argv$fq_stats, col_types = yasim_fq_stats_col_types)
fa_stats_data <- read_tsv(argv$fa_stats, col_types = yasim_fa_stats_col_types)

yasim_ground_truth <- dplyr::inner_join(fa_stats_data, fq_stats_data, by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID")) %>%
    dplyr::mutate(
        ACTUAL_RPM=ACTUAL_N_OF_READS/sum(ACTUAL_N_OF_READS)*1e6,
        ACTUAL_RPK=ACTUAL_N_OF_READS/LEN*1e3,
    ) %>%
    dplyr::mutate(
        ACTUAL_RPKM=ACTUAL_RPM/LEN*1e3,
        ACTUAL_TPM=ACTUAL_RPK/sum(ACTUAL_RPK)*1000
    )
write_tsv(yasim_ground_truth, argv$output, col_names = TRUE)

ggplot(a, aes(x=THEORETICAL_DEPTH)) +
    stat_summary(aes(y=ACTUAL_N_OF_READS), color="red") +
    stat_summary(aes(y=ACTUAL_RPM), color="blue") +
    stat_summary(aes(y=ACTUAL_RPK), color="purple") +
    stat_summary(aes(y=ACTUAL_RPKM), color="yellow") +
    stat_summary(aes(y=ACTUAL_TPM), color="green")
