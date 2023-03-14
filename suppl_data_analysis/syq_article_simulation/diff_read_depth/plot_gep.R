library("tidyverse")
library("pheatmap")

fns <- Sys.glob("ce11_as_2_isoform_depth_*.fq.stats.xz")
conditions <- fns %>%
    stringr::str_replace("ce11_as_2_isoform_depth_", "") %>%
    stringr::str_replace(".fq.stats.xz", "")

all_transcript_ids <- readr::read_tsv(
    "../ce11_as_2.gtf.gz.transcripts.tsv.xz",
    show_col_types = FALSE,
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
) %>%
    dplyr::select(TRANSCRIPT_ID)

if (file.exists("all_gep_data.parquet")) {
    all_gep_data <- arrow::read_parquet("all_gep_data.parquet")
} else {
    all_gep_data <- NULL
    for (i in seq_along(conditions)) {
        message(sprintf("%d/%d -- %s", i, length(fns), fns[i]))
        this_gep_data <- readr::read_tsv(
            fns[i],
            col_types = c(
                TRANSCRIPT_ID = col_character(),
                INPUT_DEPTH = col_double(),
                SIMULATED_N_OF_READS = col_integer(),
                TRANSCRIBED_LENGTH = col_integer(),
                SIMULATED_DEPTH = col_double()
            )
        ) %>%
            dplyr::select(TRANSCRIPT_ID, INPUT_DEPTH, SIMULATED_DEPTH) %>%
            dplyr::right_join(all_transcript_ids, by = "TRANSCRIPT_ID") %>%
            dplyr::mutate(across(where(is.numeric), replace_na, 0)) %>%
            dplyr::mutate(
                Condition = conditions[i],
                SIM_INPUT_RATIO = SIMULATED_DEPTH / INPUT_DEPTH
            )
        if (is.null(all_gep_data)) {
            all_gep_data <- this_gep_data
        } else {
            all_gep_data <- all_gep_data %>%
                dplyr::rows_append(this_gep_data)
        }
        rm(i, this_gep_data)
        gc()
    }
    arrow::write_parquet(all_gep_data, "all_gep_data.parquet")
}

g <- ggplot(all_gep_data) +
    geom_point(aes(x = log(SIM_INPUT_RATIO), y = INPUT_DEPTH), size = 0.2, alpha = 0.1) +
    facet_wrap(. ~ Condition, scales = "free") +
    theme_bw()

ggsave("gep_ratio.png", g, width = 15, height = 12)

g <- ggplot(all_gep_data) +
    geom_boxplot(aes(x=INPUT_DEPTH, y=Condition)) +
    facet_wrap(.~Condition, scales="free") +
    theme_bw()

ggsave("gep_real.png", g, width=15, height=12)

g <- ggplot(all_gep_data) +
    geom_histogram(aes(x=INPUT_DEPTH)) +
    xlim(c(0, 500)) +
    facet_wrap(.~Condition) +
    theme_bw()

ggsave("gep_real_hist.png", g, width=15, height=12)

means <- all_gep_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
        MEAN_INPUT_DEPTH = mean(INPUT_DEPTH),
        MEAN_SIMULATED_DEPTH = mean(SIMULATED_DEPTH),
        MAX_INPUT_DEPTH = max(INPUT_DEPTH),
        MAX_SIMULATED_DEPTH = mean(SIMULATED_DEPTH)
    )

print(as.data.frame(means))
