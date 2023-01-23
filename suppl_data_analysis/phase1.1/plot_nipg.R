library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

fns <- c(
    "ce11_badread_nanopore2018.fq.bam.stringtie.gtf.gene.tsv",
    "ce11_badread_nanopore2020.fq.bam.stringtie.gtf.gene.tsv",
    "ce11_badread_pacbio2016.fq.bam.stringtie.gtf.gene.tsv",
    "ce11_pbsim2_r94.fq.bam.stringtie.gtf.gene.tsv",
    "ce11_pbsim_clr.fq.bam.stringtie.gtf.gene.tsv"
)
conditions <- c(
    "NANOPORE2018",
    "NANOPORE2020",
    "PACBIO2016",
    "R94",
    "CLR"
)


#' All NIpG Data
all_nipg_data <- list()

#' Read data
for (i in seq_along(conditions)) {
    all_nipg_data[[conditions[i]]] <-
        readr::read_tsv(
            fns[i],
            show_col_types = FALSE,
            col_types = c(
                GENE_ID = col_character(),
                TRANSCRIPT_NUMBER = col_integer(),
                NAIVE_LENGTH = col_integer(),
                TRANSCRIBED_LENGTH = col_integer(),
                MAPPABLE_LENGTH = col_integer()
            )
        ) %>%
            dplyr::select(TRANSCRIPT_NUMBER) %>%
            unlist() %>%
            as.vector()
}

n_bins <- 0
for (i in conditions) {
    n_bins <- max(n_bins, max(all_nipg_data[[i]]))
}

all_nipg_data_binned <- tibble(
    x = 1:n_bins
)

for (i in conditions) {
    adding_col <- NULL
    for (n in 1:n_bins) {
        adding_col[n] <- length(which(all_nipg_data[[i]] == n))
    }
    all_nipg_data_binned <- all_nipg_data_binned %>%
        dplyr::mutate(!!i := adding_col)
}

pheatmap_data <- all_nipg_data_binned
pheatmap_data$x <- NULL

pheatmap(
    log(pheatmap_data + 1),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    main = "Logged Number of isoforms in a gene accross all conditions",
    filename = "nipg_heatmap.pdf",
    width = 8,
    height = 5
)
