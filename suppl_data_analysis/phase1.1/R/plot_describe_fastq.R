library(tidyverse)
argv <- commandArgs(trailingOnly = TRUE)

outdir_path <- argv[1]

all_tsv <- readr::read_tsv(
    file.path(outdir_path, "all.tsv"),
    col_types = cols(
        SEQID = col_character(),
        GC = col_double(),
        LEN = col_double(),
        MEANQUAL = col_double()
    ),
    progress = TRUE
)

gc_mean <- mean(all_tsv$GC)
prealign_gc_plot <- ggplot(all_tsv) +
    geom_density(aes(x = GC)) +
    xlim(c(0, 1)) +
    theme_bw() +
    ggtitle(sprintf(
        "Prealign GC plot (mean %f)",
        gc_mean
    ))
ggsave(
    file.path(outdir_path, "prealign_gc_plot.pdf"),
    prealign_gc_plot,
    height = 8,
    width = 10
)

len_mean <- mean(all_tsv$LEN)
prealign_len_plot <- ggplot(all_tsv) +
    geom_density(aes(x = LEN)) +
    theme_bw() +
    xlim(c(0, len_mean * 2)) +
    ggtitle(sprintf(
        "Prealign LEN plot (mean %f)",
        len_mean
    ))
ggsave(
    file.path(outdir_path, "prealign_len_plot.pdf"),
    prealign_len_plot,
    height = 8,
    width = 10
)

prealign_qual_plot <- ggplot(all_tsv) +
    geom_density(aes(x = MEANQUAL)) +
    theme_bw() +
    ggtitle("Prealign QUAL plot")
ggsave(
    file.path(outdir_path, "prealign_qual_plot.pdf"),
    prealign_qual_plot,
    height = 8,
    width = 10
)

extension_tsv <- readr::read_tsv(
    file.path(outdir_path, "extension_stat.tsv"),
    col_types = cols(
        POS = col_double(),
        QUAL = col_double()
    ),
    progress = TRUE
)

prealign_extension_plot <- ggplot(extension_tsv) +
    geom_line(aes(x = POS, y = QUAL), alpha=0.1) +
    xlim(c(0, NA)) +
    theme_bw() +
    ggtitle("Prealign base extension plot")
ggsave(
    file.path(outdir_path, "prealign_extension_plot.pdf"),
    prealign_extension_plot,
    height = 8,
    width = 10
)
