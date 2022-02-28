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
        SIMULATED_RPM=SIMULATED_N_OF_READS/sum(SIMULATED_N_OF_READS)*1e6,
        SIMULATED_RPK=SIMULATED_N_OF_READS/LEN*1e3,
    ) %>%
    dplyr::mutate(
        SIMULATED_RPKM=SIMULATED_RPM/LEN*1e3,
        SIMULATED_TPM=SIMULATED_RPK/sum(SIMULATED_RPK)*1000
    )
write_tsv(yasim_ground_truth, argv$output, col_names = TRUE)
