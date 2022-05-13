library("tidyverse")
library("svglite")
library("fitdistrplus")

#' Get dataset metadata
datasets_metadata <-
    readr::read_csv("datasets.csv", show_col_types = FALSE) %>%
        dplyr::mutate(filename = sprintf("./datasets/%s/%s.stringtie_guided.gtf.gene.tsv", Model, Read))

#' Raw sequence run accessions
all_run_accessions <- unlist(datasets_metadata$Read)

#' All GEP data in depth
all_gep_data <- list()

#' Read data
for (i in seq_len(nrow(datasets_metadata))) {
    data_tuple <- datasets_metadata[i,]
    all_gep_data[[data_tuple$Read]] <-
        readr::read_tsv(data_tuple$filename, show_col_types = FALSE) %>%
            dplyr::select(TRANSCRIPT_NUMBER) %>%
            unlist() %>%
            as.vector()
}

n_bins <- 0
for (i in all_run_accessions) {
    n_bins <- max(n_bins, max(all_gep_data[[i]]))
}

all_gep_data_reformatted <- tibble(
    x= 1:n_bins
)

for (i in all_run_accessions) {
    adding_col <- NULL
    for (n in 1:n_bins){
        adding_col[n] <- length(which(all_gep_data[[i]]==n))
    }
    all_gep_data_reformatted <- all_gep_data_reformatted %>%
        dplyr::mutate(!!i := adding_col)
}

all_gep_data_reformatted_long <- all_gep_data_reformatted %>%
    tidyr::gather(key = "Read", value = "GeneNumber", -x)

g <- ggplot(all_gep_data_reformatted_long) +
    geom_line(aes(x = x, y=GeneNumber, color = Read), stat="identity") +
    scale_x_discrete(name = "N. Isoforms") +
    scale_y_continuous(name = "log(N. Genes)", trans="log") +
    theme_bw() +
    ggtitle("Number of isoforms in a gene accross all samples")
ggsave("gep_n_isoforms.svg", g, width = 8, height = 5)

all_gep_data_reformatted_mean <- all_gep_data_reformatted %>%
    dplyr::mutate(count=as.integer(rowMeans(all_gep_data_reformatted)))

fit_common <- function(data) {
    fnbiom <- fitdistrplus::fitdist(as.integer(data$count), "nbinom", method = "mme")
    fpois <- fitdistrplus::fitdist(as.integer(data$count), "pois", method = "mme")
    fgamma <- fitdistrplus::fitdist(data$count, "gamma", method = "mme")
    flnorm <- fitdistrplus::fitdist(data$count + 1, "lnorm")
    fexp <- fitdistrplus::fitdist(data$count, "exp", method = "mme")
    # fweibull <- fitdistrplus::fitdist(data$count, "weibull", method = "mle")
    len_data <- dim(data)[1]
    data <- data %>%
        dplyr::mutate(
            lnorm_theo = sort(rlnorm(
                len_data,
                meanlog = flnorm$estimate[1],
                sdlog = flnorm$estimate[2]
            )) - 1,
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
            ))# ,
            # weibull_theo = sort(rweibull(
            #     len_data,
            #     scale = fweibull$estimate[1],
            #     shape = fweibull$estimate[2]
            # ))
        )
    return(
        list(
            fnbiom = fnbiom,
            fpois = fpois,
            fgamma = fgamma,
            flnorm = flnorm,
            fexp = fexp,
            # fweibull = fweibull,
            data = data
        )
    )
}
plot_density <- function(retl) {
    data <- retl$data
    g <- ggplot(data) +
        geom_density(aes(x = count), color = "red", size = 3) +
        geom_density(aes(x = lnorm_theo), color = "blue") +
        geom_density(aes(x = gamma_theo), color = "purple") +
        geom_density(aes(x = nbiom_theo), color = "black") +
        geom_density(aes(x = pois_theo), color = "green") +
        geom_density(aes(x = exp_theo), color = "yellow") +
        scale_x_continuous(
            limits = c(1E-3, exp(n_bins)),
            trans="log"
        ) +
        scale_y_continuous(
            limits = c(0, 0.3)
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
        scale_x_continuous(limits = c(0, n_bins + 1)) +
        scale_y_continuous(limits = c(0, n_bins + 1)) +
        theme_bw()
}
gs_all <- fit_common(all_gep_data_reformatted_mean)
g1 <- plot_density(gs_all)
ggsave("fitted_n_isoforms.svg",g1, width=8, height=5)
write_csv(all_gep_data_reformatted_mean, "fitted_n_isoforms.csv")

#' Get number of records
n_of_records <- nrow(all_gep_data)
