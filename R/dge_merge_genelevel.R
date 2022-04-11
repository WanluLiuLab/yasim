file_description <- "dge_merge_genelevel.R -- Merge TSVs Produced by Gene-Level Quantifiers for Downstream DGE Analysis"

#' This file is used to merge TSVs produced by Transcript-Level Quantifiers for Downstream DGE Analysis.
#' TODO: Gene level unmapped

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(
    p, "--libfile", type = "character",
    help = "Directory of libfiles to read."
)
p <- add_argument(
    p, "--featureCounts_tsv", default = c(), nargs = Inf,
    help = ".tsv file produced by featureCounts"
)
p <- add_argument(
    p, "--salmon_quant_sf", default = c(), nargs = Inf,
    help = ".sf file produced by salmon quant."
)
p <- add_argument(
    p, "--cpptetgs_tsv", default = c(), nargs = Inf,
    help = "*.tsv produced by CPPTETGS"
)
p <- add_argument(
    p, "--htseq_count_tsv", default = c(), nargs = Inf,
    help = "*.tsv produced by htseq-count"
)
p <- add_argument(
    p, "--fq_stats", default = c(), nargs = Inf,
    help = "*.fq.stats produced by YASIM simulator"
)
p <- add_argument(
    p, "--unmapped_stats", default = c(), nargs = Inf,
    help = "*_unmapped.stats produced by aligner"
)
p <- add_argument(
    p, "--fa_stats", type = "character",
    help = "*.fa.stats produced by YASIM transcriptor"
)
p <- add_argument(p, "--output", help = "Output basename", type = "character")
argv <- parse_args(p)

source(file.path(argv$libfile, "loader.R"))
source(file.path(argv$libfile, "cpptetgs.R"))
source(file.path(argv$libfile, "featureCounts.R"))
source(file.path(argv$libfile, "htseq.R"))
source(file.path(argv$libfile, "salmon.R"))
source(file.path(argv$libfile, "yasim.R"))

#' The metadata of transcripts, like length, position, GC, etc.
fa_stats_data <- read_tsv(argv$fa_stats, col_types = yasim_fa_stats_col_types)

#' The entire large wide table
all_table <- fa_stats_data %>%
    dplyr::group_by(GENE_ID) %>%
    dplyr::select(GENE_ID) %>%
    dplyr::summarize(LEN = 1)


#' Number of reads in ground truth.
sum_actual_n_of_reads <- c()

#' Number of mapped reads in ground truth.
sum_simulated_n_of_reads <- c()

#' Function that calculates TPM, RPM, RPK and RPKM from Number of Reads and Transcript Length
#' @param all_table: Table that need to be parsed. Should be some dataframe-like structure.
#' @param step_name: Name of the step. Can be name of Simulators (Ground Truth, in this case) or Quantifiers.
#' @param n: The replication number.
mutate_tables_for_rpkm <- function(all_table, step_name, n) {
    all_table <- all_table %>%
        mutate(across(where(is.numeric), replace_na, 0)) %>%
        dplyr::filter(LEN != 0) # Filter transcripts added by accident.
    all_table <- all_table %>%
        dplyr::mutate(
            !!sprintf("%s_%d_ACTUAL_LOG2", step_name, n) :=
                log2(.[[!!sprintf("%s_%d_ACTUAL_N_OF_READS", step_name, n)]] + 1)
        )

    return(all_table)
}

#' Start parsing ground truth.

if (!is.na(argv$fq_stats)) {
    n <- 1
    for (fq_stats in argv$fq_stats) {
        fq_stats_data <- aggregate_fq_stats_data_into_gene_level(
            get_fq_stats_data(fq_stats, n),
            fa_stats_data,
            n
        )
        all_table <- dplyr::full_join(
            all_table,
            fq_stats_data,
            by = "GENE_ID"
        )
        sum_simulated_n_of_reads <- c(
            sum_simulated_n_of_reads,
            sum(all_table[[sprintf("%s_%d_ACTUAL_N_OF_READS", "YASIM", n)]])
        )
        message(sprintf("Loaded %d reads from %s", sum_simulated_n_of_reads[n], fq_stats))
        all_table <- mutate_tables_for_rpkm(all_table, "YASIM", n)
        n <- n + 1
    }
}
# if (!is.na(argv$unmapped_stats)) {
#     n <- 1
#     for (unmapped_stats in argv$unmapped_stats) {
#         unmapped_stats_data <- get_unmapped_stats_data(unmapped_stats, n)
#         all_table <- dplyr::full_join(all_table, unmapped_stats_data, by = TRANSCRIPT_ID)
#         all_table <- replace(all_table, is.na(all_table), 0)
#         sum_actual_n_of_reads <- c(
#             sum_actual_n_of_reads,
#             sum(all_table[[sprintf("%s_%d_ACTUAL_N_OF_READS", "YASIM_MAPPED", n)]])
#         )
#         message(sprintf("Loaded %d reads from %s", sum_actual_n_of_reads[n], fq_stats))
#         all_table <- mutate_tables_for_rpkm(all_table, "YASIM_MAPPED", n)
#         n <- n + 1
#     }
# }
if (!is.na(argv$featureCounts_tsv)) {
    n <- 1
    for (featureCounts_tsv in argv$featureCounts_tsv) {
        featureCounts_data <- get_featureCounts_data(
            featureCounts_tsv,
            n,
            "GENE_ID"
        )

        all_table <- dplyr::full_join(
            all_table,
            featureCounts_data,
            by = "GENE_ID"
        )

        all_table <- mutate_tables_for_rpkm(all_table, "FEATURECOUNTS", n)
        n <- n + 1
    }
}
if (!is.na(argv$salmon_quant_sf)) {
    n <- 1
    for (salmon_quant_sf in argv$salmon_quant_sf) {
        salmon_quant_sf_data <- get_salmon_data(
            salmon_quant_sf,
            n,
            "GENE_ID"
        )

        all_table <- dplyr::full_join(
            all_table,
            salmon_quant_sf_data,
            by = "GENE_ID"
        )

        all_table <- mutate_tables_for_rpkm(all_table, "SALMON", n)
        n <- n + 1
    }
}
if (!is.na(argv$htseq_count_tsv)) {
    n <- 1
    for (htseq_count_tsv in argv$htseq_count_tsv) {
        htseq_count_tsv_data <- get_htseq_count_data(
            htseq_count_tsv,
            n,
            "GENE_ID"
        )

        all_table <- dplyr::full_join(
            all_table,
            htseq_count_tsv_data,
            by = "GENE_ID"
        )

        all_table <- mutate_tables_for_rpkm(all_table, "HTSEQ_COUNT", n)
        n <- n + 1
    }
}
if (!is.na(argv$cpptetgs_tsv)) {
    n <- 1
    for (cpptetgs_tsv in argv$cpptetgs_tsv) {
        cpptetgs_tsv_data <- get_cpptetgs_data(
            cpptetgs_tsv,
            n,
            "GENE_ID"
        )

        all_table <- dplyr::full_join(
            all_table,
            cpptetgs_tsv_data,
            by = "GENE_ID"
        )

        all_table <- mutate_tables_for_rpkm(all_table, "CPPTETGS", n)
        n <- n + 1
    }
}
write_tsv(all_table, paste0(argv$output, ".tsv"))
