library("tidyverse")
library("pheatmap")
library("ggridges")

actual_fns <- Sys.glob("ce11_as_2_isoform_depth_*.fq.stats")
conditions <- actual_fns %>%
    stringr::str_replace("ce11_as_2_isoform_depth_", "") %>%
    stringr::str_replace(".fq.stats", "")


gene_id_transcript_id <- readr::read_tsv(
    "ce11_as_2.gtf.transcripts.tsv",
    show_col_types = FALSE,
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
) %>%
    dplyr::select(TRANSCRIPT_ID, GENE_ID)
all_transcript_ids <- gene_id_transcript_id %>%
    dplyr::select(TRANSCRIPT_ID)


if (file.exists("all_gep_data.parquet")) {
    all_gep_data <- arrow::read_parquet("all_gep_data.parquet")
} else {
    all_gep_data <- NULL
    for (i in seq_along(conditions)) {
        message(sprintf("%d/%d -- %s", i, length(actual_fns), actual_fns[i]))
        this_gep_data <- readr::read_tsv(
            actual_fns[i],
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
    geom_point(aes(y = log10(SIMULATED_DEPTH), x = log10(INPUT_DEPTH)), size = 0.2, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    ylab("Log10 LLRG Simulated") +
    xlab("Log10 LLRG Input") +
    ggtitle("LLRG Error Summarize") +
    facet_wrap(. ~ Condition) +
    theme_bw()

ggsave("gep_ratio.png", g, width = 15, height = 12)

g <- ggplot(all_gep_data) +
    geom_histogram(aes(x = SIMULATED_DEPTH)) +
    xlim(c(0, 500)) +
    ylab("N. Events") +
    xlab("LLRG Simulated Depth") +
    ggtitle("LLRG Simulated Depth Histogram") +
    facet_wrap(. ~ Condition, scales = "free") +
    theme_bw()

ggsave("gep_real_hist.pdf", g, width = 15, height = 12)

means <- all_gep_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
        MEAN_INPUT_DEPTH = mean(INPUT_DEPTH),
        MEAN_SIMULATED_DEPTH = mean(SIMULATED_DEPTH),
        MAX_INPUT_DEPTH = max(INPUT_DEPTH),
        MAX_SIMULATED_DEPTH = mean(SIMULATED_DEPTH)
    )


g <- ggplot(means) +
    geom_bar(aes(x = MEAN_SIMULATED_DEPTH, y = Condition), stat = "identity") +
    ylab("Mean LLRG Simulated Depth") +
    xlab("Condition") +
    ggtitle("Mean Simulated Depth bar Plot") +
    theme_bw()

ggsave("gep_real_bar.pdf", g, width = 15, height = 12)

g <- ggplot(all_gep_data) +
    geom_density_ridges_gradient(aes(
        x = SIMULATED_DEPTH + 1,
        y = Condition
    )) +
    ylab("density") +
    xlim(c(0, 200)) +
    theme_ridges()

ggsave("gep_real_ridges.pdf", g, width = 10, height = 8)

all_gep_data_with_gene_id <- all_gep_data %>%
    dplyr::inner_join(gene_id_transcript_id, by = "TRANSCRIPT_ID") %>%
    dplyr::select(GENE_ID, TRANSCRIPT_ID, SIMULATED_DEPTH, Condition)

if (file.exists("gep_inside_gene_isoform_level_variation.parquet")) {
    gep_inside_gene_isoform_level_variation <- arrow::read_parquet(
        "gep_inside_gene_isoform_level_variation.parquet"
    )
} else {
    gep_inside_gene_isoform_level_variation <- data.frame()
    for (gene_id in sample(unique(all_gep_data_with_gene_id$GENE_ID), 1000)) {
        for (condition in unique(all_gep_data_with_gene_id$Condition)) {
            this_gep_data_with_gene_id <- all_gep_data_with_gene_id %>%
                dplyr::filter(GENE_ID == gene_id, Condition == condition)
            this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
                dplyr::cross_join(this_gep_data_with_gene_id)
            gep_inside_gene_isoform_level_variation <- rbind(
                gep_inside_gene_isoform_level_variation,
                data.frame(
                    var = this_gep_data_with_gene_id$SIMULATED_DEPTH.x / this_gep_data_with_gene_id$SIMULATED_DEPTH.y,
                    Condition = condition
                )
            )
        }
    }
    arrow::write_parquet(
        gep_inside_gene_isoform_level_variation,
        "gep_inside_gene_isoform_level_variation.parquet"
    )
}

g <- ggplot(gep_inside_gene_isoform_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    xlab("Variance (Log 10 Fold Change)") +
    ggtitle("Distribution of Variation Inside a Gene") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()
ggsave("gep_var.pdf", g, width = 15, height = 12)

if (file.exists("gep_isoform_level_variation.parquet")) {
    gep_isoform_level_variation <- arrow::read_parquet("gep_isoform_level_variation.parquet")
} else {
    gep_isoform_level_variation <- data.frame()
    for (condition in unique(all_gep_data_with_gene_id$Condition)) {
        this_gep_data_with_gene_id <- all_gep_data_with_gene_id %>%
            dplyr::filter(Condition == condition) %>%
            head(3000)
        this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
            dplyr::cross_join(this_gep_data_with_gene_id)
        gep_isoform_level_variation <- rbind(
            gep_isoform_level_variation,
            data.frame(
                var = this_gep_data_with_gene_id$SIMULATED_DEPTH.x / this_gep_data_with_gene_id$SIMULATED_DEPTH.y,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_isoform_level_variation, "gep_isoform_level_variation.parquet")
}

g <- ggplot(gep_isoform_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    xlab("Variance (Log 10 Fold Change)") +
    ggtitle("Distribution of Variation among all Isoforms") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave("gep_isoform_level_variation.pdf", g, width = 15, height = 12)

if (file.exists("gep_gene_level_variation.parquet")) {
    gep_gene_level_variation <- arrow::read_parquet("gep_gene_level_variation.parquet")
} else {
    gep_inside_gene_isoform_level_variation3 <- data.frame()
    for (condition in unique(all_gep_data_with_gene_id$Condition)) {
        print(condition)
        this_gep_data_with_gene_id <- all_gep_data_with_gene_id %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::group_by(GENE_ID) %>%
            dplyr::summarise(SIMULATED_DEPTH = mean(SIMULATED_DEPTH)) %>%
            unique() %>%
            head(n = 1000)
        this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
            dplyr::cross_join(this_gep_data_with_gene_id)
        gep_gene_level_variation <- rbind(
            gep_gene_level_variation,
            data.frame(
                var = this_gep_data_with_gene_id$SIMULATED_DEPTH.x / this_gep_data_with_gene_id$SIMULATED_DEPTH.y,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_gene_level_variation, "gep_gene_level_variation.parquet")
}

g <- ggplot(gep_gene_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    xlab("Variance (Log 10 Fold Change)") +
    ggtitle("Distribution of Variation among all Genes") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave("gep_gene_level_variation.pdf", g, width = 15, height = 12)
