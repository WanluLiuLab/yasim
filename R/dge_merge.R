file_description <- ""

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(
    p, "--featureCounts_tsv", default = c(), nargs = Inf,
    help = ".tsv file produced by featureCounts"
)
p <- add_argument(
    p, "--salmon_quant_sf", default = c(), nargs = Inf,
    help = ".sf file produced by salmon quant."
)
p <- add_argument(
    p, "--stringtie_quant_tsv", default = c(), nargs = Inf,
    help = ".tsv file produced by parse_stringtie_into_tsv.py"
)
p <- add_argument(
    p, "--cpptetgs_tsv", default = c(), nargs = Inf,
    help = "*.tsv produced by CPPTETGS"
)
p <- add_argument(
    p, "--fq_stats", default = c(), nargs = Inf,
    help = "*.fq.stats produced by YASIM simulator"
)
p <- add_argument(p, "--fa_stats", help = "*.fa.stats produced by YASIM transcriptor", type = "character")
p <- add_argument(p, "--output", help = "Output basename", type = "character")
argv <- parse_args(p)

source(argv$libfile)

fa_stats_data <- read_tsv(argv$fa_stats, col_types = yasim_fa_stats_col_types)
#' The entire large wide table
all_table <- fa_stats_data

sum_actual_n_of_reads <- c()

mutate_tables_for_rpkm <- function(all_table, step_name, n, sum_actual_n_of_reads) {
    all_table <- all_table %>%
        dplyr::mutate(
            !!sprintf("%s_%d_ACTUAL_RPM", step_name, n) :=
                .[[!!sprintf("%s_%d_ACTUAL_N_OF_READS", step_name, n)]] /
                    sum_actual_n_of_reads[n] * 1e6,
            !!sprintf("%s_%d_ACTUAL_RPK", step_name, n) :=
                .[[!!sprintf("%s_%d_ACTUAL_N_OF_READS", step_name, n)]] / LEN * 1e3,
        )
    # print(names(all_table))
    sum_rpk <- sum(
        replace_na(all_table[[sprintf("%s_%d_ACTUAL_RPK", step_name, n)]], 0)
    )
    # write_tsv(all_table,paste0(argv$output, ".tsv"))
    all_table <- all_table %>% dplyr::mutate(
        !!sprintf("%s_%d_ACTUAL_RPKM", step_name, n) :=
            .[[!!sprintf("%s_%d_ACTUAL_RPM", step_name, n)]] / LEN * 1e3,
        !!sprintf("%s_%d_ACTUAL_TPM", step_name, n) :=
            .[[!!sprintf("%s_%d_ACTUAL_RPK", step_name, n)]] / sum_rpk * 1e3
    )
    return(all_table)
}

if (!is.na(argv$fq_stats)) {
    n <- 1
    for (fq_stats in argv$fq_stats) {
        fq_stats_data <- get_yasim_data(fq_stats, n)
        all_table <- dplyr::full_join(all_table, fq_stats_data, by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID"))
        all_table <- replace(all_table, is.na(all_table), 0)
        sum_actual_n_of_reads <- c(
            sum_actual_n_of_reads,
            sum(all_table[[sprintf("%s_%d_ACTUAL_N_OF_READS", "YASIM", .GlobalEnv$n)]])
        )
        message(sprintf("Loaded %d reads from %s", sum_actual_n_of_reads[n], fq_stats))
        all_table <- mutate_tables_for_rpkm(all_table, "YASIM", n, sum_actual_n_of_reads)
        n <- n + 1
    }
}
if (!is.na(argv$featureCounts_tsv)) {
    n <- 1
    for (featureCounts_tsv in argv$featureCounts_tsv) {
        featureCounts_data <- get_featureCounts_data(featureCounts_tsv, n)

        all_table <- dplyr::full_join(
            all_table,
            featureCounts_data,
            by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID")
        )
        all_table <- mutate(all_table, across(where(is.numeric), replace_na, 0))
        all_table <- mutate_tables_for_rpkm(all_table, "FEATURECOUNTS", n, sum_actual_n_of_reads)
        n <- n + 1
    }
}
if (!is.na(argv$salmon_quant_sf)) {
    n <- 1
    for (salmon_quant_sf in argv$salmon_quant_sf) {
        salmon_quant_sf_data <- get_salmon_data(salmon_quant_sf, n)

        all_table <- dplyr::full_join(
            all_table,
            salmon_quant_sf_data,
            by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID")
        )
        all_table <- mutate(all_table, across(where(is.numeric), replace_na, 0))
        all_table <- mutate_tables_for_rpkm(all_table, "SALMON", n, sum_actual_n_of_reads)
        n <- n + 1
    }
}
if (!is.na(argv$stringtie_quant_tsv)) {
    n <- 1
    for (stringtie_quant_tsv in argv$stringtie_quant_tsv) {
        stringtie_quant_tsv_data <- get_stringtie_data(stringtie_quant_tsv, n)

        all_table <- dplyr::full_join(
            all_table,
            stringtie_quant_tsv_data,
            by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID")
        )
        all_table <- mutate(all_table, across(where(is.numeric), replace_na, 0))
        # No mutate, done in reader.
        n <- n + 1
    }
}
if (!is.na(argv$cpptetgs_tsv)) {
    n <- 1
    for (cpptetgs_tsv in argv$cpptetgs_tsv) {
        cpptetgs_tsv_data <- get_cpptetgs_data(cpptetgs_tsv, n)

        all_table <- dplyr::full_join(
            all_table,
            cpptetgs_tsv_data,
            by = c("TRANSCRIPT_ID" = "TRANSCRIPT_ID")
        )
        all_table <- mutate(all_table, across(where(is.numeric), replace_na, 0))
        all_table <- mutate_tables_for_rpkm(all_table, "CPPTETGS", n, sum_actual_n_of_reads)
        n <- n + 1
    }
}
write_tsv(all_table, paste0(argv$output, ".tsv"))
