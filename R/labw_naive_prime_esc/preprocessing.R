#!/usr/bin/env Rscript
#' Extracted from DETE.R
rm(list = ls())
is_debug <- FALSE
file_description <- "preprocessing.R -- Pre-processing step of DETE"


#' Function to load a package while supressing all their outputs.
#' A simple wrapper for library()
#' @param name The name of the package to use
load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

# Rprof("Profile.out")
message(file_description)
message(cat(commandArgs(), sep = ' '))
message("================ Loading packages ================")
load_package("tibble")
load_package("svglite")
load_package("argparser")
load_package("tidyverse")
load_package("ggpubr")

# si <- sessionInfo()
# pid <- Sys.getpid()
# save.image(si,file=paste0(pid, "_sessionInfo.RData"))

lazy_loader <- function(final_filename) {
    message(sprintf("Lazy loading from %s", final_filename))
    return(read_csv(final_filename,
                    col_names = TRUE,
                    col_types = cols(
                        name = col_character(),
                        FINE2a = col_double(),
                        FINE2b = col_double(),
                        TesR7A = col_double(),
                        TesR7B = col_double())
    ))
}

merge_fc_data <- function(FINE2a_filename,
                          FINE2b_filename,
                          TesR7A_filename,
                          TesR7B_filename,
                          suffix) {
    final_filename <- paste0(suffix, ".csv")
    if (file.exists(final_filename)) {
        matrix_raw <- lazy_loader(final_filename)
    } else {
        col_spec <- cols(
            Geneid = col_character(),
            Chr = col_character(),
            Start = col_integer(),
            End = col_integer(),
            Strand = col_character(),
            Length = col_integer(),
            count = col_double()
        )
        col_name <- c(
            "Geneid",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
            "count"
        )
        FINE2a_f <- read_tsv(FINE2a_filename, col_names = col_name, comment = "#", col_types = col_spec)
        FINE2b_f <- read_tsv(FINE2b_filename, col_names = col_name, comment = "#", col_types = col_spec)
        TesR7A_f <- read_tsv(TesR7A_filename, col_names = col_name, comment = "#", col_types = col_spec)
        TesR7B_f <- read_tsv(TesR7B_filename, col_names = col_name, comment = "#", col_types = col_spec)
        message("Read finished, merging...")
        matrix_raw <- cbind(FINE2a_f[, c(1, 7)], FINE2b_f[, 7], TesR7A_f[, 7], TesR7B_f[, 7])
        colnames(matrix_raw) <- c("name", "FINE2a", "FINE2b", "TesR7A", "TesR7B")
        matrix_raw <- matrix_raw %>%
            dplyr::filter(FINE2a + FINE2b + TesR7A + TesR7B != 0)
        write_csv(matrix_raw, final_filename)
    }
    return(matrix_raw)
}

merge_tetgs_data <- function(FINE2a_filename,
                             FINE2b_filename,
                             TesR7A_filename,
                             TesR7B_filename,
                             suffix) {
    final_filename <- paste0(suffix, ".csv")
    if (file.exists(final_filename)) {
        matrix_raw <- lazy_loader(final_filename)
    } else {
        col_spec <- cols(
            TE_NAME = col_character(),
            count = col_double()
        )
        col_name <- c("TE_NAME", "count")
        FINE2a_raw <- read_tsv(FINE2a_filename, col_names = col_name, comment = "#", col_types = col_spec)
        FINE2b_raw <- read_tsv(FINE2b_filename, col_names = col_name, comment = "#", col_types = col_spec)
        TesR7A_raw <- read_tsv(TesR7A_filename, col_names = col_name, comment = "#", col_types = col_spec)
        TesR7B_raw <- read_tsv(TesR7B_filename, col_names = col_name, comment = "#", col_types = col_spec)
        message("Read finished, merging...")
        matrix_raw <- dplyr::full_join(FINE2a_raw, FINE2b_raw, by = "TE_NAME")
        matrix_raw <- dplyr::full_join(matrix_raw, TesR7A_raw, by = "TE_NAME")
        matrix_raw <- dplyr::full_join(matrix_raw, TesR7B_raw, by = "TE_NAME")
        colnames(matrix_raw) <- c("name", "FINE2a", "FINE2b", "TesR7A", "TesR7B")
        matrix_raw <- replace(matrix_raw, is.na(matrix_raw), 0)
        write_csv(matrix_raw, final_filename)
    }
    return(matrix_raw)
}

message("================ Loading data ================")
if (!is_debug) {
    p <- arg_parser(file_description)
    p <- add_argument(p, "--wd", help = "Working directory", type = "character", default = "..")
    p <- add_argument(p, "--FINE2a", help = "FINE2a", type = "character")
    p <- add_argument(p, "--FINE2b", help = "FINE2b", type = "character")
    p <- add_argument(p, "--TesR7A", help = "TesR7A", type = "character")
    p <- add_argument(p, "--TesR7B", help = "TesR7B", type = "character")
    p <- add_argument(p, "--suffix", help = "suffix", type = "character")
    p <- add_argument(p, "--is_fc", help = "Whether the data is from FC", flag = TRUE)
    argv <- parse_args(p)
    setwd(argv$wd)
    suffix <- argv$suffix
    if (argv$is_fc) {
        readcount <- merge_fc_data(argv$FINE2a, argv$FINE2b, argv$TesR7A, argv$TesR7B, suffix)
    } else {
        readcount <- merge_tetgs_data(argv$FINE2a, argv$FINE2b, argv$TesR7A, argv$TesR7B, suffix)
    }
} else {
    message("Debug mode")
    suffix <- "fc"
    readcount <- lazy_loader("fc.csv")
}

message("================ Quality Control ================")

#' Get number of genes/isoforms covered in each sample
gene_covered_df <- readcount %>%
    dplyr::summarise_if(is.numeric, function(x) sum(x > 0)) %>%
    t() %>%
    as_tibble(rownames = "name", .name_repair = function(names) { "value" })

gene_covered_plot <- ggbarplot(
    gene_covered_df, x = "name", y = "value",
    title = sprintf("Total genes covered (%s)", suffix),
    xlab = "Sample Name",
    ylab = "No. of Gene Covered"
)

ggsave(sprintf("%s_total_gene_covered.tiff", suffix),
       plot = gene_covered_plot)

#' Get copies of genes/isoforms covered in each sample
total_coverage_df <- readcount %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    t() %>%
    as_tibble(rownames = "name", .name_repair = function(names) { "value" })
total_coverage_plot <- ggbarplot(
    total_coverage_df, x = "name", y = "value",
    title = sprintf("log10(total counts over genes (%s))", suffix),
    xlab = "Sample Name",
    ylab = "Coverage of Sample"
)
ggsave(sprintf("%s_total_coverage.tiff", suffix),
       plot = total_coverage_plot)
