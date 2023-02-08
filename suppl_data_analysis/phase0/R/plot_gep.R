library("tidyverse")
library("ggpubr")
library("gamlss")
library("fitdistrplus")
library("parallel")
library("rmutil")
library("ggpmisc")

cl <- parallel::makeCluster(parallel::detectCores() - 1)
clusterExport(cl, varlist = ls())

all_data <- read_tsv(
    "../all_data.tsv",
    show_col_types = FALSE,
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        SEQNAME = col_character(),
        START = col_double(),
        END = col_double(),
        STRAND = col_character(),
        ABSOLUTE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        GC = col_double(),
        PACB_AVG_DEPTH = col_double(),
        NANOPORE_AVG_DEPTH = col_double(),
        ILLM_AVG_DEPTH = col_double()
    )
)

scale_depth <- function(vec) {
    (vec - min(vec)) / (max(vec) - min(vec)) * 5000 + 1
}

all_data_nano <- all_data %>%
    dplyr::rename(count = NANOPORE_AVG_DEPTH) %>%
    dplyr::filter(count > 0)
all_data_pacb <- all_data %>%
    dplyr::rename(count = PACB_AVG_DEPTH) %>%
    dplyr::filter(count > 0)
all_data_illm <- all_data %>%
    dplyr::rename(count = ILLM_AVG_DEPTH) %>%
    dplyr::filter(count > 0)

all_data_nano$count <- scale_depth(all_data_nano$count)
all_data_pacb$count <- scale_depth(all_data_pacb$count)
all_data_illm$count <- scale_depth(all_data_illm$count)

mutate_data <- function(data) {
    dplyr::mutate(data,
                  rank = dplyr::dense_rank(dplyr::desc(count)),
                  rank_zipfs = rank^(-2),
                  zipfs_freq = max(count) / rank,
                  rank_l10n = log10(rank),
                  rank_zipfs_l10n = log10(rank_zipfs),
                  zipfs_freq_l10n = log10(zipfs_freq),
                  count_l10n = log10(count),
    )
}

plot_zipf <- function(data) {
    ggplot(data) +
        geom_point(aes(x = count_l10n, y = rank_zipfs_l10n), color = "blue") +
        geom_point(aes(x = zipfs_freq_l10n, y = rank_zipfs_l10n), color = "red") +
        geom_smooth(aes(x = count_l10n, y = rank_zipfs_l10n), method = 'lm') +
        theme_bw()
}

all_data_nano <- mutate_data(all_data_nano)
g <- plot_zipf(all_data_nano) +
    ggtitle(
        "Zipf's Distribution Fitting for Nanopore Data ",
        "Blue dots: Actual data. Red dots: Theoritical Distribution. Line: Fitted Distribution"
    )
ggsave("gep_zipf_nanopore.pdf", g, width = 10, height = 8)

all_data_pacb <- mutate_data(all_data_pacb)
g <- plot_zipf(all_data_pacb) +
    ggtitle(
        "Zipf's Distribution Fitting for PacBio Data ",
        "Blue dots: Actual data. Red dots: Theoritical Distribution. Line: Fitted Distribution"
    )
ggsave("gep_zipf_pacbio.pdf", g, width = 10, height = 8)

all_data_illm <- mutate_data(all_data_illm)
g <- plot_zipf(all_data_illm) +
    ggtitle(
        "Zipf's Distribution Fitting for Illumina Data ",
        "Blue dots: Actual data. Red dots: Theoritical Distribution. Line: Fitted Distribution"
    )
ggsave("gep_zipf_illumina.pdf", g, width = 10, height = 8)

fit_common <- function(data) {
    fnbiom <- fitdistrplus::fitdist(as.integer(data$count), "nbinom", method = "mme")
    fpois <- fitdistrplus::fitdist(as.integer(data$count), "pois", method = "mme")
    fgamma <- fitdistrplus::fitdist(data$count, "gamma", method = "mme")
    flnorm <- fitdistrplus::fitdist(data$count, "lnorm", method = "mme")
    fexp <- fitdistrplus::fitdist(data$count, "exp", method = "mme")
    fweibull <- fitdistrplus::fitdist(data$count, "weibull", method = "mle")
    len_data <- dim(data)[1]
    data <- data %>%
        dplyr::arrange(desc(rank))
    data <- data %>%
        dplyr::mutate(
            lnorm_theo = sort(rlnorm(
                len_data,
                meanlog = flnorm$estimate[1],
                sdlog = flnorm$estimate[2]
            )),
            pois_theo = sort(rpois(
                len_data,
                lambda = fpois$estimate[1]
            )),
            gamma_theo = sort(rgamma(
                len_data,
                shape = fgamma$estimate[1],
                rate = fgamma$estimate[2]
            )),
            nbiom_theo = sort(rnbinom(
                len_data,
                size = fnbiom$estimate[1],
                mu = fnbiom$estimate[2]
            )),
            exp_theo = sort(rexp(
                len_data,
                rate = fexp$estimate[1]
            )),
            weibull_theo = sort(rweibull(
                len_data,
                scale = fweibull$estimate[1],
                shape = fweibull$estimate[2]
            ))
        )
    return(
        list(
            fnbiom = fnbiom,
            fpois = fpois,
            fgamma = fgamma,
            flnorm = flnorm,
            fexp = fexp,
            fweibull = fweibull,
            data = data
        )
    )
}

plot_density <- function(retl) {
    data <- retl$data

    g <- ggplot(data) +
        geom_density(aes(x = log10(count)), color = "red", size = 3) +
        geom_density(aes(x = log10(lnorm_theo)), color = "blue") +
        geom_density(aes(x = log10(gamma_theo)), color = "purple") +
        geom_density(aes(x = log10(nbiom_theo)), color = "black") +
        scale_x_continuous(
            limits = c(-1, NA)
        )
    gb <- ggplot_build(g)
    g +
        stat_peaks(
            data = gb[['data']][[1]],
            aes(x = x, y = density),
            geom = "text"
        ) +
        theme_bw()

}

plot_qq <- function(retl) {
    data <- retl$data
    ggplot(data, aes(sample = count)) +
        geom_qq(distribution = qlnorm, dparams = retl$flnorm$estimate, color = "blue") +
        geom_qq(distribution = qgamma, dparams = retl$fgamma$estimate, color = "purple") +
        geom_qq(distribution = qnbinom, dparams = retl$fnbiom$estimate, color = "black") +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw()
}

gs_nano <- fit_common(all_data_nano)
g <- plot_density(gs_nano) +
    ggtitle(
        "Fitting of Commonly Seen Distributions for Nanopore Data",
        "Blue: LogNormal. Purple: Gamma. Black: Negative Biomial. Red: Real"
    )
ggsave("gep_common_nanopore.pdf", g, width = 10, height = 8)

g <- plot_qq(gs_nano) +
    ggtitle(
        "Fitting of Commonly Seen Distributions for Nanopore Data",
        "Blue: LogNormal. Purple: Gamma. Black: Negative Biomial. Red: Real"
    )
ggsave("gep_common_nanopore_qq.pdf", g, width = 10, height = 8)
gs_pacb <- fit_common(all_data_pacb)
plot_density(gs_pacb)
plot_qq(gs_pacb)
gs_illm <- fit_common(all_data_illm)
plot_density(gs_illm)
plot_qq(gs_illm)
