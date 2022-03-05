library("tidyverse")
library("ggpubr")
library("gamlss")
library("fitdistrplus")

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
        NANOPORE_AVG_DEPTH>10,
        PACB_AVG_DEPTH>10,
        NANOPORE_AVG_DEPTH+PACB_AVG_DEPTH != Inf
    )

fig_data_depth <- all_data %>%
     tidyr::gather(
         key="key",
         value = "value",
         -TRANSCRIPT_ID,
         -GENE_ID,
         -SEQNAME,
         -START,
         -END,
         -STRAND,
         -LEN,
         -GC
     )

ggplot(fig_data_depth, aes(x=log(value))) + geom_histogram(aes(color=key))

diff_data <- all_data %>%
    dplyr::mutate(DIFF=(NANOPORE_AVG_DEPTH-PACB_AVG_DEPTH)/(NANOPORE_AVG_DEPTH+PACB_AVG_DEPTH)) %>%
    dplyr::filter(abs(DIFF)<0.5)

ggplot(diff_data, aes(x=PACB_AVG_DEPTH, y=NANOPORE_AVG_DEPTH)) + geom_point()


fit <- fitDist(
    diff_data$NANOPORE_AVG_DEPTH,
    k = 2,
    type = "realAll",
    try.gamlss = TRUE,
    parallel="snow"
)

plot(fit)
summary(fit)
