library("tidyverse")
library("pheatmap")
library("ggridges")

fns <- Sys.glob("ce11_as_*.gtf.gene.tsv")
conditions <- fns %>%
    stringr::str_replace("ce11_as_", "") %>%
    stringr::str_replace(".gtf.gene.tsv", "")

fns <- c(
    "../ce11.ncbiRefSeq.gtf.gene.tsv",
    "../ce11_as_2.gtf.gene.tsv",
    fns
)

conditions <- c(
    "REFERENCE",
    "2",
    conditions
)

if (file.exists("all_nipg_data_binned.parquet")) {
    all_nipg_data_binned <- arrow::read_parquet("all_nipg_data_binned.parquet")
    all_nipg_data_long <- arrow::read_parquet("all_nipg_data_long.parquet")
} else {
    all_nipg_data <- list()
    all_nipg_data_long <- data.frame()
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
            all_nipg_data_long <- rbind(
                all_nipg_data_long,
                data.frame(
                    TRANSCRIPT_NUMBER=all_nipg_data[[conditions[i]]],
                    Condition=conditions[i]
                )
            )
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
    arrow::write_parquet(all_nipg_data_long, "all_nipg_data_long.parquet")
}

pheatmap_data <- all_nipg_data_binned
pheatmap_data$x <- NULL

pheatmap(
    log(pheatmap_data + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Logged Number of Isoforms in a Gene accross all Conditions",
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

g <- ggplot(all_nipg_data_long, aes(x = TRANSCRIPT_NUMBER, y = Condition)) +
    geom_boxplot() +
    stat_summary(fun.x="mean",color="red", shape=13) +
    scale_x_continuous("Number of Isoforms in a Gene", limits = c(0, 20)) +
    theme_bw()
ggsave("nipg_hist2.pdf", g, width = 8, height = 5)

# g <- ggplot(all_nipg_data_binned_long) +
#     geom_bar(aes(x = x, y = Count), stat = "identity") +
#     theme_bw() +
#     facet_grid(Condition~., scale = "free_y")
# ggsave("nipg_hist3.pdf", g, width = 8, height = 5)

dir.create("nipg_hist4")

for (n in unique(all_nipg_data_binned_long$Condition)){
    this_nipg_data_binned_long <- all_nipg_data_binned_long %>%
        dplyr::filter(Condition==n)
    this_all_nipg_data_long <- all_nipg_data_long %>%
        dplyr::filter(Condition==n)
    mean_bar <- mean(this_all_nipg_data_long$TRANSCRIPT_NUMBER)
    g <- ggplot(this_nipg_data_binned_long) +
        geom_bar(aes(x = x, y = Count), stat = "identity") +
        geom_vline(xintercept=mean_bar, color="red") +
        scale_x_continuous("Number of Isoforms in a Gene", limits = c(0, 20)) +
        theme_bw()
    ggsave(file.path("nipg_hist4", sprintf("%s.pdf", n)), g, width = 15, height = 3)
}

all_nipg_data_long %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(TRANSCRIPT_NUMBER=mean(TRANSCRIPT_NUMBER)) %>%
    print()

# g <- ggplot(all_nipg_data_long, aes(x = TRANSCRIPT_NUMBER, y = Condition)) +
#     geom_point(position="jitter", alpha=0.2, size=0.2) +
#     geom_violin(color="red") +
#     scale_x_continuous("Number of Isoforms Number in a Gene", limits = c(0, 20)) +
#     theme_bw()
# ggsave("nipg_violin.pdf", g, width = 8, height = 5)

#' Count number of genes and transcripts in each sample

fns <- Sys.glob("ce11_as_*.gtf.transcripts.tsv")

conditions <- fns %>%
    stringr::str_replace("ce11_as_", "") %>%
    stringr::str_replace(".gtf.transcripts.tsv", "")

fns <- c(
    "../ce11.ncbiRefSeq.gtf.transcripts.tsv",
    "../ce11_as_2.gtf.transcripts.tsv",
    fns
)

conditions <- c(
    "REFERENCE",
    "2",
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

g1 <- ggplot(all_nipg_ngenes_nisoforms) +
    geom_bar(aes(y = Condition, x = nGenes), stat = "identity") +
    scale_x_continuous("Total Number of Genes") +
    scale_y_discrete("Targeted Complexity") +
    theme_bw()

g2 <- ggplot(all_nipg_ngenes_nisoforms) +
    geom_bar(aes(y = Condition, x = nIsoforms), stat = "identity") +
    scale_x_continuous("Total Number of Isoforms") +
    scale_y_discrete("") +
    theme_bw()

g <- cowplot::plot_grid(g1, g2, nrow = 1, labels = LETTERS[1:4]) +
    ggtitle("Total Number of Genes and Isoforms in each Simulation")

ggsave("nipg_ngenes_nisoforms.pdf", g, width = 8, height = 5)
