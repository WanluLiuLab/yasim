library("tidyverse")
library("pheatmap")

fns <- Sys.glob("ce11_as_*.gtf.gz.gene.tsv.xz")
conditions <- fns %>%
    stringr::str_replace("ce11_as_", "") %>%
    stringr::str_replace(".gtf.gz.gene.tsv.xz", "")

fns <- c(
    "../ce11.ncbiRefSeq.gtf.gz.gene.tsv.xz",
    "../ce11_as_2.gtf.gz.gene.tsv.xz",
    fns
)

conditions <- c(
    "REFERENCE",
    "2",
    conditions
)

if (file.exists("all_nipg_data_binned.parquet")) {
    all_nipg_data_binned <- arrow::read_parquet("all_nipg_data_binned.parquet")
} else {
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
    arrow::write_parquet(all_nipg_data_binned, "all_nipg_data_binned.parquet")
}

pheatmap_data <- all_nipg_data_binned
pheatmap_data$x <- NULL

pheatmap(
    log(pheatmap_data + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Logged Number of isoforms in a gene accross all conditions",
    filename = "nipg_heatmap.pdf",
    width = 8,
    height = 5
)

all_nipg_data_binned_long <- all_nipg_data_binned %>%
    tidyr::gather(key = "Condition", value = "Count", -x)

g <- ggplot(all_nipg_data_binned_long) +
    geom_bar(aes(x = x, y = Count), stat = "identity") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    facet_wrap(. ~ Condition)
ggsave("nipg_hist.pdf", g, width = 8, height = 5)

#' Count number of genes and transcripts in each sample

fns <- Sys.glob("ce11_as_*.gtf.gz.transcripts.tsv.xz")
conditions <- fns %>%
    stringr::str_replace("ce11_as_", "") %>%
    stringr::str_replace(".gtf.gz.transcripts.tsv.xz", "")

fns <- c(
    "ce11.ncbiRefSeq.gtf.gz.transcripts.tsv.xz",
    fns
)

conditions <- c(
    "REFERENCE",
    conditions
)


if (file.exists("all_nipg_ngenes_nisoforms.parquet")) {
    all_nipg_ngenes_nisoforms <- arrow::read_parquet("all_nipg_ngenes_nisoforms.parquet")
} else {
        #' All NIpG Data
    all_data <- NULL

        #' Read data
    for (i in seq_along(conditions)) {
        this_data <- readr::read_tsv(
            fns[i],
            show_col_types = FALSE,
            col_types = c(
                TRANSCRIPT_ID = col_character(),
                GENE_ID = col_character(),
                TRANSCRIPT_NUMBER = col_integer(),
                NAIVE_LENGTH = col_integer(),
                TRANSCRIBED_LENGTH = col_integer(),
                EXON_NUMBER = col_integer()
            )
        ) %>%
            dplyr::select(TRANSCRIPT_ID, GENE_ID) %>%
            dplyr::mutate(Condition = conditions[i])
        if (is.null(all_data)) {
            all_data <- this_data
        } else {
            all_data <- all_data %>%
                dplyr::rows_append(this_data)
        }
    }
    all_nipg_ngenes_nisoforms <- all_data %>%
        dplyr::group_by(Condition) %>%
        dplyr::transmute(
            nGenes = length(unique(GENE_ID)),
            nIsoforms = length(unique(TRANSCRIPT_ID)),
            Condition = Condition
        ) %>%
        dplyr::distinct()

    arrow::write_parquet(all_nipg_ngenes_nisoforms, "all_nipg_ngenes_nisoforms.parquet")
}

g <- ggplot(all_nipg_ngenes_nisoforms) +
    geom_bar(aes(y = Condition, x = nGenes), stat = "identity") +
    theme_bw()

ggsave("nipg_ngenes.pdf", g, width = 8, height = 5)

g <- ggplot(all_nipg_ngenes_nisoforms) +
    geom_bar(aes(y = Condition, x = nIsoforms), stat = "identity") +
    theme_bw()

ggsave("nipg_nisoforms.pdf", g, width = 8, height = 5)
