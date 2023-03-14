library("tidyverse")
library("ggridges")
library("arrow")

all_data <- NULL

fns <- Sys.glob("ce11_*pbsim3*.maf.gz.rlen.tsv")
conditions <- fns %>%
    stringr::str_replace("ce11_", "") %>%
    stringr::str_replace(".maf.gz.rlen.tsv", "")

if (file.exists("all_last_aligned_maf_data.parquet")){
    all_data <- arrow::read_parquet("all_last_aligned_maf_data.parquet")
} else{
    
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
    message(sprintf("Processing %d/%d", i, length(conditions)))
    rm(this_data, i)
    gc()
}
arrow::write_parquet(all_data, "all_last_aligned_maf_data.parquet")

}
simulated_transcript_stats <- readr::read_tsv(
    "../ce11_as_2.gtf.gz.transcripts.tsv.xz",
    show_col_types = FALSE
)
all_data_sim <- all_data %>%
    dplyr::inner_join(
        simulated_transcript_stats,
        by = c("ALIGNED_TRANSCRIPT_ID" = "TRANSCRIPT_ID")
    ) %>%
    dplyr::sample_n(10000) %>%
    dplyr::transmute(
        Condition = Condition,
        TRANSCRIPT_ID = SIMULATED_TRANSCRIPT_ID,
        READ_COMPLETENESS = READ_LENGTH / TRANSCRIBED_LENGTH
    )

g <- ggplot(all_data_sim) +
    geom_density_ridges(
        aes(
            x = READ_COMPLETENESS,
            y = Condition
        ),
        alpha = 0
    ) +
    ylab("density") +
    xlim(c(0.2, 1.5)) +
    theme_ridges() +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness.pdf", g, width = 8, height = 8)

g <- ggplot(all_data_sim) +
    geom_boxplot(
        aes(
            x = READ_COMPLETENESS,
            y = Condition
        ),
        alpha = 0
    ) +
    ylab("density") +
    theme_ridges() +
    ggtitle("Read Completeness of all conditions")

ggsave("last_read_completeness_box.pdf", g, width = 8, height = 8)
