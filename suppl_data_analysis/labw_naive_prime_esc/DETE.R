#!/usr/bin/env Rscript
#' Originally written by SU Yaqi (DETE_ORI.R)
#' Reformatted by YU Zhejian
rm(list = ls())
is_debug <- FALSE
file_description <- "DETE.R -- Analyse differentially expressed TEs"


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
load_package("BiocParallel")
load_package("tibble")
load_package("DESeq2")
load_package("pheatmap")
load_package("argparser")
load_package("tidyverse")
load_package("ggpubr")

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

message("================ Loading data ================")
if (!is_debug) {
    p <- arg_parser(file_description)
    p <- add_argument(p, "--wd", help = "Working directory", type = "character", default = "..")
    p <- add_argument(p, "--suffix", help = "suffix", type = "character")
    argv <- parse_args(p)
    setwd(argv$wd)
    suffix <- argv$suffix
    readcount <- lazy_loader(sprintf("%s.csv", suffix))
} else {
    message("Debug mode")
    suffix <- "fc"
    readcount <- lazy_loader("fc.csv")
}

message("================ Loading Data into DESeq2 ================")
#' Transfer the data type into integer
#' DESeq2 accept DataFrame (S4) or data.frame (S3) with pure integer only,
#' so need to transform the tibble into dataframe
#' and remove rownames
readcount_int <- readcount %>%
    dplyr::mutate_if(is.numeric, as.integer) %>%
    as.data.frame()
rownames(readcount_int) <- readcount_int$name
readcount_int$name <- NULL
configure <- tibble(condition = factor(c("FINE", "FINE", "TesR", "TesR")),
                    type = c("r1", "r2", "r1", "r2"))
dds <- DESeqDataSetFromMatrix(countData = readcount_int,
                              colData = configure,
                              design = ~condition)

message("================ Normalization ================")
#' Perform VST, rlog and log2(x)+1 normalization
#' TODO: Select between VST and rlog
vst_normalized_model <- vst(dds, blind = FALSE)
vst_normalized_data <- assay(vst_normalized_model)
rlog_normalized_model <- rlog(dds, blind = FALSE)
rlog_normalized_data <- assay(rlog_normalized_model)
dds <- estimateSizeFactors(dds)
#' TODO: Why? is it a RPKM?
#' divide the counts by the size factors or normalization factors before returning
log2p1_normalized_data <- log2(counts(dds, normalized = TRUE) + 1)

#' Plot normalization plot
normalization_df <- dplyr::bind_rows(
    as_tibble(log2p1_normalized_data[, 1:2]) %>%
        dplyr::mutate(normalization_type = "log2(x + 1)"),
    as_tibble(vst_normalized_data[, 1:2]) %>%
        dplyr::mutate(normalization_type = "vst"),
    as_tibble(rlog_normalized_data[, 1:2]) %>%
        dplyr::mutate(normalization_type = "rlog")
)
colnames(normalization_df)[1:2] <- c("x", "y")
normalization_plot <- ggplot(normalization_df, aes(x = x, y = y)) +
    geom_hex(bins = 80) +
    coord_fixed() +
    facet_grid(. ~ normalization_type) +
    labs(xlab = "X",
         ylab = "Y",
         title = sprintf("Normalization method (%s)", suffix)
    )
ggsave(sprintf("%s_normalized.tiff", suffix),
       plot = normalization_plot)

message("================ PCA ================")
#' Generate the PCA plots
pca_plot <- plotPCA(vst_normalized_model)
ggsave(sprintf("%s_pca.tiff", suffix),
       plot = pca_plot)

message("================ Calling DETE ================")
#' Pre filtering
#' While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful:
#' reduce the memory size of the dds data object
#' increase the speed of the transformation and testing functions within DESeq2.
#' Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.
#' Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.
TE_to_keep <- rowSums(counts(dds)) >= 10
dds <- dds[TE_to_keep,]

#' Call DETE
#' The use of parallel in DESeq() or result() is deprecated due to the reason that
#' parallelization makes the original function slower.
parallel_cluster <- MulticoreParam()
dds <- DESeq(dds)

base_mean_DETE <- counts(dds, normalized = TRUE) %>%
    as_tibble(rownames = "name") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        baseMean_FINE = mean(c(FINE2a, FINE2b)),
        baseMean_TESR = mean(c(TesR7A, TesR7B))
    ) %>%
    dplyr::select(name, baseMean_FINE, baseMean_TESR)

write_csv(base_mean_DETE, sprintf("%s_base.csv", suffix))
#' Set the p value and log2FC threshold
DETE_noFiltered <- as_tibble(results(dds, lfcThreshold = 2, alpha = 0.05), rownames = "name")
write_csv(DETE_noFiltered, sprintf("%s_DETE_noFiltered.csv", suffix))

#' MA plot
tiff(sprintf("%s_ma.tiff", suffix))
plotMA(dds)
dev.off()

message("================ Saving Significant DETEs ================")
#' Get the DEGs with padj<0.05 and |log2FoldChange| >=2 as output.
DETE_filtered <- DETE_noFiltered %>%
    dplyr::filter(
        abs(log2FoldChange) >= 2,
        padj < 0.05
    )
write_csv(DETE_filtered, sprintf("%s_DETE.csv", suffix))

pheatmap_data <- vst_normalized_data %>%
    as.data.frame() %>%
    dplyr::mutate(var = stats::var(c(FINE2a, FINE2b, TesR7A, TesR7B))) %>%
    dplyr::arrange(desc(var)) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::select(!var)
heatmap_plot <- pheatmap(pheatmap_data)
ggsave(sprintf("%s_heatmap.tiff", suffix),
       plot = heatmap_plot)
pheatmap_data <- pheatmap_data %>% as_tibble(rownames = "name")
write_csv(pheatmap_data, sprintf("%s_DETE_top20.csv", suffix))
# unlink(paste0(pid, "_sessionInfo.RData"))
# Rprof(NULL)
