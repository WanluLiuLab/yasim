library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")

all_data <- NULL

fns <- Sys.glob("ce11_*.fq_trans.maf.rlen.tsv")
conditions <- fns %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".fq_trans.maf.rlen.tsv", "")

for (i in seq_along(fns)) {
    this_data <- readr::read_tsv(
        fns[i],
        col_types = c(
            ALIGNED_TRANSCRIPT_ID = col_character(),
            SIMULATED_TRANSCRIPT_ID = col_character(),
            READ_LENGTH = col_integer()
        ),
        progress = TRUE,
        quote = "\'"
    ) %>%
        dplyr::mutate(
            Condition = conditions[i]
        )
    if (is.null(all_data)) {
        all_data <- this_data
    } else {
        all_data <- all_data %>%
            dplyr::rows_append(this_data)
    }
    rm(this_data, i)
}
aligned_transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.chr1.gtf.transcripts.tsv",
    show_col_types = FALSE
)
simulated_transcript_stats <- readr::read_tsv(
    "ce11.ncbiRefSeq.chr1_as.gtf.transcripts.tsv",
    show_col_types = FALSE
)
correctness <- all_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
        correctness = sum(SIMULATED_TRANSCRIPT_ID == ALIGNED_TRANSCRIPT_ID) / n()
    ) %>%
    dplyr::ungroup()

g <- ggplot(correctness) +
    geom_bar(aes(y = Condition, x = correctness), stat = "identity") +
    xlim(c(0.0, 1.0)) +
    theme_bw()
ggsave("last_correctness.pdf", g, width = 5, height = 4)

all_data_sim <- all_data %>%
    dplyr::inner_join(
        simulated_transcript_stats,
        by = c("SIMULATED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::transmute(
        Condition = sprintf("%s_sim", Condition),
        TRANSCRIPT_ID = SIMULATED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )
all_data_aln <- all_data %>%
    dplyr::inner_join(
        aligned_transcript_stats,
        by = c("ALIGNED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::transmute(
        Condition = sprintf("%s_aln", Condition),
        TRANSCRIPT_ID = ALIGNED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )

all_data_long <- dplyr::rows_append(all_data_sim, all_data_aln)

g <- ggplot(all_data_long) +
    geom_density_ridges_gradient(
        aes(
            x = READ_COMPLETENESS,
            y = Condition
        )
    ) +
    ylab("density") +
    xlim(c(0.7, 1)) +
    theme_ridges() +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness.pdf", g, width = 8, height = 8)

all_data_public <- all_data %>%
    dplyr::inner_join(
        aligned_transcript_stats,
        by = c("ALIGNED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::transmute(
        ALIGNED_READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH,
        Condition = Condition,
        READ_LENGTH = READ_LENGTH,
        SIMULATED_TRANSCRIPT_ID = SIMULATED_TRANSCRIPT_ID
    ) %>%
    dplyr::inner_join(
        simulated_transcript_stats,
        by = c("SIMULATED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::transmute(
        SIMULATED_READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH,
        Condition = Condition,
        ALIGNED_READ_COMPLETENESS = ALIGNED_READ_COMPLETENESS
    )

all_data_p <- NULL

all_conditions <- unique(all_data_public$Condition)

for (accession in all_conditions){
    this_data <- all_data_public %>%
                    dplyr::filter(Condition==accession)
    distance <- dist(
            rbind(
                density(this_data$SIMULATED_READ_COMPLETENESS, from = 0, to = 1)$y,
                density(this_data$ALIGNED_READ_COMPLETENESS, from = 0, to = 1)$y
            ),
            method = "euclidean"
        )[1]
    pvalue <- wilcox.test(
        density(this_data$SIMULATED_READ_COMPLETENESS, from = 0, to = 1)$y,
        density(this_data$ALIGNED_READ_COMPLETENESS, from = 0, to = 1)$y
    )$p.value
    this_data_p <- tibble(
        Condition=accession,
        pvalue = pvalue,
        distance =distance
    )
    if (is.null(all_data_p)){
        all_data_p <- this_data_p
    } else{
        all_data_p <- dplyr::rows_append(
            all_data_p, this_data_p
        )
    }
    rm(accession, this_data)
}

ggplot(all_data_p) +
    geom_bar(aes(y=Condition, x=pvalue), stat="identity") +
    theme_bw()
