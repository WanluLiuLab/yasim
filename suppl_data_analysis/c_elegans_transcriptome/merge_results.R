library(tidyverse)
library(ggpubr)
library(gamlss)
library(fitdistrplus)
library(parallel)

setwd("suppl_data_analysis/c_elegans_transcriptome")

depth_data_col_type <- cols(
    TRANSCRIPT_ID=col_character(),
    BASE=col_number(),
    DEPTH=col_number()
)

depth_data_col_name <- c(
    "TRANSCRIPT_ID",
    "BASE",
    "DEPTH"
)

yasim_fa_stats_col_types <- cols(
    TRANSCRIPT_ID = col_character(),
    GENE_ID = col_character(),
    SEQNAME = col_character(),
    START = col_number(),
    END = col_number(),
    STRAND = col_character(),
    LEN = col_number(),
    GC = col_number()
)

nanopore_depth_data <- read_tsv(
    "nanopore_transcriptome_1.bam.depth.tsv",
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
)%>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(NANOPORE_AVG_DEPTH=mean(DEPTH))
pbsim_depth_data <- read_tsv(
    "pacbio_transcriptome_1.bam.depth.tsv",
    col_types = depth_data_col_type,
    col_names = depth_data_col_name
)%>%
    dplyr::group_by(TRANSCRIPT_ID) %>%
    dplyr::summarise(PACB_AVG_DEPTH=mean(DEPTH))

fa_stats_data <- read_tsv(
    "ce11.reference_transcripts.fa.stats",
    col_types = yasim_fa_stats_col_types
)

all_data <- fa_stats_data %>%
    dplyr::full_join(pbsim_depth_data, by="TRANSCRIPT_ID") %>%
    dplyr::full_join(nanopore_depth_data, by="TRANSCRIPT_ID") %>%
    dplyr::mutate(across(where(is.numeric), replace_na, 0)) %>%
    dplyr::filter(
        NANOPORE_AVG_DEPTH+PACB_AVG_DEPTH != Inf
    )

write_tsv(all_data, "all_data.tsv")


NANOPORE_AVG_DEPTH <- all_data$NANOPORE_AVG_DEPTH %>% .[! . == 0]
PACB_AVG_DEPTH <- all_data$PACB_AVG_DEPTH %>% .[! . == 0]

fit_gamma_ont <- fitdistrplus::fitdist(NANOPORE_AVG_DEPTH, "gamma")
fit_nb_ont <- fitdistrplus::fitdist(as.integer(NANOPORE_AVG_DEPTH), "nbinom")
fit_gamma_pacb <- fitdistrplus::fitdist(PACB_AVG_DEPTH, "gamma")
fit_nb_pacb <- fitdistrplus::fitdist(as.integer(PACB_AVG_DEPTH), "nbinom")

fit_gamma_ont$aic > fit_nb_ont$aic
fit_gamma_pacb$aic > fit_nb_pacb$aic

fit_gamma_ont$bic > fit_nb_ont$bic
fit_gamma_pacb$bic > fit_nb_pacb$bic


cl <- parallel::makeCluster(parallel::detectCores()-1)

clusterExport(cl, varlist = ls())

psr <- parSapply(cl, 1:1e6, function(x){
    c(
        fitdistrplus::fitdist(as.integer(
            sample(NANOPORE_AVG_DEPTH, as.integer(length(NANOPORE_AVG_DEPTH)/10))
        ), "nbinom")$estimate,
        fitdistrplus::fitdist(as.integer(
            sample(PACB_AVG_DEPTH, as.integer(length(PACB_AVG_DEPTH)/10))
        ), "nbinom")$estimate
    )
})


