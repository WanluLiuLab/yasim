#!/usr/bin/env Rscript
#' Originally written by XIANG Xinyu
#' Adopted by SU Yaqi (Plots_code_ORI.R)
#' Reformatted by YU Zhejian
#'
#' TODO: The 80 and 8 should be parameterized
rm(list = ls())
file_description <- "Plots_codes.R -- Create additional plots for DETE.R"
is_debug <- FALSE

load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

# Rprof("Profile.out")
message(file_description)
message(cat(commandArgs(), sep = ' '))
message("================ Loading packages ================")
load_package("tibble")
load_package("gridExtra")
load_package("grid")
load_package("tidyverse")
load_package("ggpubr")
load_package("argparser")


if (!is_debug) {
    p <- arg_parser(file_description)
    p <- add_argument(p, "--wd", help = "Working directory", type = "character", default = ".")
    p <- add_argument(p, "--suffix", help = "suffix", type = "character")
    argv <- parse_args(p)
    setwd(argv$wd)
    suffix <- argv$suffix
} else {
    suffix <- "fc"
}


message("================ Reading TE references ... ================")
TE_family_count <- read_tsv(
    "TE_family_count.tsv",
    col_types = cols(
        gene_id = col_character(),
        total_count = col_integer()
    ),
    col_names = TRUE
)
TE_gene_transcript_class <- read_tsv(
    'TE_gene_transcript_class.tsv',
    col_types = cols(
        gene_id = col_character(),
        transcript_id = col_character(),
        class_id = col_character()
    ),
    col_names = TRUE
)
DETE_noFiltered <- read_csv(
    sprintf("%s_DETE_noFiltered.csv", suffix),
    col_types = cols(
        name = col_character(),
        baseMean = col_double(),
        log2FoldChange = col_double(),
        lfcSE = col_double(),
        stat = col_double(),
        pvalue = col_double(),
        padj = col_double()
    ),
    col_names = TRUE
)
base_mean_DETE <- read_csv(
    sprintf("%s_base.csv", suffix),
    col_types = cols(
        name = col_character(),
        baseMean_FINE = col_double(),
        baseMean_TESR = col_double()
    ),
    col_names = TRUE
)
readcount <- read_csv(
    sprintf("%s.csv", suffix), ,
    col_names = TRUE,
    col_types = cols(
        name = col_character(),
        FINE2a = col_double(),
        FINE2b = col_double(),
        TesR7A = col_double(),
        TesR7B = col_double())
)
message("================ Joining Data ... ================")
DETE_noFiltered <- dplyr::inner_join(DETE_noFiltered, TE_gene_transcript_class, by = c("name" = "transcript_id"))
base_mean_DETE <- dplyr::inner_join(base_mean_DETE, TE_gene_transcript_class, by = c("name" = "transcript_id"))
#' DETE with fold change >= 2
DETE_up_foldchange2 <- dplyr::filter(DETE_noFiltered, log2FoldChange >= 1)
#' DETE with fold change <= 0.5
DETE_dw_foldchange2 <- dplyr::filter(DETE_noFiltered, log2FoldChange <= -1)


message("================ Producing Pie Chart ... ================")
#' The portion of identified DETE, up- and down-regulation
n_DETE_up_foldchange2 <- nrow(DETE_up_foldchange2)
n_DETE_dw_foldchange2 <- nrow(DETE_dw_foldchange2)
sum_DETE_foldchange2 <- n_DETE_up_foldchange2 + n_DETE_dw_foldchange2

pie_chart_labels <- c(
    sprintf("Up-regulated TE copies\nn=%d, %.2f%%", n_DETE_up_foldchange2, 100 * n_DETE_up_foldchange2 / sum_DETE_foldchange2),
    sprintf("Down-regulated TE copies\nn=%d, %.2f%%", n_DETE_dw_foldchange2, 100 * n_DETE_dw_foldchange2 / sum_DETE_foldchange2)
)
pie_chart_data <- tibble(value = c(n_DETE_up_foldchange2, n_DETE_dw_foldchange2))
pie_chart_plot <- ggpie(pie_chart_data,
                        x = "value",
                        label = pie_chart_labels,
                        lab.pos = "in",
                        fill = c("light blue", "light yellow"))

ggsave(sprintf("%s_pie_chart.svg", suffix),
       plot = pie_chart_plot,
       device = "svg")

message("================ Producing Bar Plot ... ================")
#' The portion of identified isoforms in all possible isoforms
DETE_up_foldchange2_gene_id_data <- dplyr::count(DETE_up_foldchange2, gene_id)
DETE_dw_foldchange2_gene_id_data <- dplyr::count(DETE_dw_foldchange2, gene_id)
#' Only kept TE gene_id with more than 8 isoforms of up/down regulated are loaded.
#' Only TE with more than 80 copies are loaded
#' Ranked by the difference of portion of identified isoforms in all possible isoforms
OEvsWT_gene_id_data <- dplyr::full_join(
    DETE_up_foldchange2_gene_id_data,
    DETE_dw_foldchange2_gene_id_data,
    by = "gene_id", suffix = c("up", "dw")) %>%
    dplyr::rename(up = nup, dw = ndw) %>%
    tidyr::replace_na(list(up = 0.0, dw = 0.0)) %>%
    dplyr::inner_join(TE_family_count, by = "gene_id") %>%
            #' TODO: No idea why have this filter
    dplyr::filter(total_count >= 80) %>%
    dplyr::mutate(
        up_percent = up / total_count,
        dw_percent = dw / total_count
    ) %>%
    dplyr::mutate(diff = up_percent - dw_percent)

OEvsWT_gene_id_data_top10 <- OEvsWT_gene_id_data %>%
            #' TODO: No isea why have this filter
    dplyr::filter(up >= 8 | dw >= 8) %>%
    dplyr::arrange(desc(diff)) %>%
    dplyr::filter(between(row_number(), 1, 10) | between(row_number(), n() - 9, n()))
# print(OEvsWT_gene_id_data_top10)
#' Factorize for plotting
OEvsWT_gene_id_data_factored <- OEvsWT_gene_id_data_top10 %>%
    dplyr::mutate(gene_id = factor(gene_id, levels = rev(as.character(gene_id))))

graph_colors <- colorRampPalette(c("dodgerblue4", 'indianred4'))(20)
middle_axis_plot <- ggplot(OEvsWT_gene_id_data_factored, aes(x = 1, y = gene_id)) +
    geom_text(aes(label = gene_id)) +
    ggtitle("") +
    ylab(NULL) +
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(color = NA),
          axis.ticks.x = element_line(color = NA),
          plot.margin = unit(c(1, -1, 1, -1), "mm"))

subplot_theme <- theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_line(colour = NA),
          panel.grid.major = element_line(colour = NA),
          plot.margin = unit(c(1, -1, 1, 0), "mm"))

up_plot <- ggplot(data = OEvsWT_gene_id_data_factored, aes(x = gene_id, y = up_percent, fill = gene_id)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = graph_colors) +
    ggtitle(sprintf("Up regulated TEs in TesR (vs FINE) %s", suffix)) +
    subplot_theme +
    scale_y_reverse(limits = c(0.3, 0), breaks = seq(0.3, 0, -0.1)) +
    coord_flip()

dw_plot <- ggplot(data = OEvsWT_gene_id_data_factored, aes(x = gene_id, y = dw_percent, fill = gene_id)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = graph_colors) +
    ggtitle(sprintf("Down regulated TEs in TesR (vs FINE) %s", suffix)) +
    subplot_theme +
    scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    coord_flip()

arranged_plot <- grid.arrange(
    ggplot_gtable(ggplot_build(up_plot)),
    ggplot_gtable(ggplot_build(middle_axis_plot)),
    ggplot_gtable(ggplot_build(dw_plot)),
    ncol = 3, widths = c(4 / 9, 1 / 9, 4 / 9))

ggsave(sprintf("%s_top10_bar_plot.svg", suffix),
       plot = arranged_plot,
       device = "svg")

message("================ Producing Top 3 regulated Scatter Plot ... ================")
#' Plot isoforms of top 3 difference of portion of identified isoforms in all possible isoforms
#' Isoforms of abs(log2foldchange) >= 1 is colored.
#' axis: log2(baseMean + 1)
DETE_top3_gene_id <- dplyr::slice(OEvsWT_gene_id_data_top10, 1:3, 18:20)$gene_id
DETE_log2fc_transcript_id <- dplyr::filter(DETE_noFiltered, abs(log2FoldChange) >= 1)$name

base_mean_DETE_top3 <- base_mean_DETE %>%
    dplyr::filter(gene_id %in% DETE_top3_gene_id) %>%
    dplyr::mutate(
        Type = dplyr::if_else(
            name %in% DETE_log2fc_transcript_id,
            gene_id,
            "Nochange"
        )
    ) %>%
    dplyr::mutate(
        Type = factor(Type, levels = c(as.character(DETE_top3_gene_id), "Nochange")),
        baseMean_FINE = log2(baseMean_FINE + 1),
        baseMean_TESR = log2(baseMean_TESR + 1)
    )

middle_line <- geom_abline(slope = 1,
                           intercept = 0,
                           colour = "black",
                           linetype = "dotdash")
scatter_theme <- theme_bw() +
    theme(axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 15),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA))
top3_scatter_plot <- ggplot(base_mean_DETE_top3) +
    geom_point(data = dplyr::filter(base_mean_DETE_top3, Type == "Nochange"),
               mapping = aes(y = baseMean_FINE, x = baseMean_TESR),
               color = "grey") +
    geom_point(data = dplyr::filter(base_mean_DETE_top3, Type != "Nochange"),
               mapping = aes(y = baseMean_FINE, x = baseMean_TESR, colour = Type)) +
    middle_line +
    scatter_theme +
    labs(title = sprintf("DETE_noFiltered (%s)", suffix),
         y = "log2(basemean+1) in FINE",
         x = "log2(basemean+1) in TesR") +
    scale_colour_manual(values = c(
        "#490A3D",
        "#BD1550",
        "#E97F02",
        "#F8CA00",
        "#8A9B0F",
        "#3878AA")) +
    xlim(0, 8) +
    ylim(0, 8)

ggsave(sprintf("%s_top3_scatter_plot.svg", suffix),
       plot = top3_scatter_plot,
       device = "svg",
       width = 9, height = 8
)

message("================ Producing Expression Level Scatter Plot ... ================")
#' size: The portion of up/down regulated isoforms that
#' axis: log2(CPM+1)

DETE_cpm_mean <- readcount %>%
            #' TODO: No idea what this line do
    dplyr::filter(name %in% DETE_noFiltered$name) %>%
    dplyr::inner_join(TE_gene_transcript_class, by = c("name" = "transcript_id")) %>%
            #' Calculate isoform-level CPM (Count per Million)
    dplyr::mutate(
        FINE2a_cpm = FINE2a / sum(FINE2a) * 10^6,
        FINE2b_cpm = FINE2b / sum(FINE2b) * 10^6,
        TesR7A_cpm = TesR7A / sum(TesR7A) * 10^6,
        TesR7B_cpm = TesR7B / sum(TesR7B) * 10^6
    ) %>%
            #' Calculate gene-level CPM
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(
        FINE_cpm_mean = mean(c(FINE2a_cpm, FINE2b_cpm)),
        TesR_cpm_mean = mean(c(TesR7A_cpm, TesR7B_cpm))
    ) %>%
    dplyr::inner_join(OEvsWT_gene_id_data, by = "gene_id") %>%
    dplyr::arrange(desc(diff)) %>%
    dplyr::mutate(
        FINE_cpm_mean = log2(FINE_cpm_mean + 1),
        TesR_cpm_mean = log2(TesR_cpm_mean + 1),
        top5 = "noSig",
        top5 = dplyr::if_else(between(row_number(), 1, 5), "upTop5", top5),
        top5 = dplyr::if_else(between(row_number(), n() - 4, n()), "dwTop5", top5)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        size = max(c(up_percent, dw_percent)),
        up_percent = NULL,
        dw_percent = NULL
    )
scatter_plot <- ggplot(DETE_cpm_mean, aes(FINE_cpm_mean, TesR_cpm_mean)) +
    scatter_theme +
    labs(title = sprintf("DETE in %s", suffix),
         y = "log2(CPM+1) in FINE",
         x = "log2(CPM+1) in TesR") +
    geom_point(aes(size = size, colour = top5), alpha = 0.4) +
    geom_text(data = dplyr::filter(DETE_cpm_mean, top5 != "noSig"),
              aes(label = gene_id), size = 1.5) +
    middle_line +
    xlim(0, NA) +
    ylim(0, NA)
ggsave(sprintf("%s_scatter_plot.svg", suffix),
       plot = scatter_plot,
       device = "svg",
       width = 7, height = 6
)

