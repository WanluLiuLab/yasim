library("tidyverse")
library("pheatmap")
library("ggridges")

actual_fns <- Sys.glob("ce11_as_2_isoform_depth_*.tsv")
conditions <- actual_fns %>%
    stringr::str_replace("ce11_as_2_isoform_depth_", "") %>%
    stringr::str_replace(".tsv", "")


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


if (file.exists("all_gep_data_simulated.parquet")) {
    all_gep_data_simulated <- arrow::read_parquet("all_gep_data_simulated.parquet")
} else {
    all_gep_data_simulated <- NULL
    for (i in seq_along(conditions)) {
        message(sprintf("%d/%d -- %s", i, length(actual_fns), actual_fns[i]))
        this_gep_data <- readr::read_tsv(
            actual_fns[i],
            col_types = c(
                TRANSCRIPT_ID = col_character(),
                DEPTH = col_double()
            )
        ) %>%
            dplyr::select(TRANSCRIPT_ID, DEPTH) %>%
            dplyr::right_join(all_transcript_ids, by = "TRANSCRIPT_ID") %>%
            dplyr::mutate(across(where(is.numeric), \(x) replace_na(x, 0))) %>%
            dplyr::mutate(Condition=conditions[i])
        if (is.null(all_gep_data_simulated)) {
            all_gep_data_simulated <- this_gep_data
        } else {
            all_gep_data_simulated <- all_gep_data_simulated %>%
                dplyr::rows_append(this_gep_data)
        }
        rm(i, this_gep_data)
        gc()
    }
    arrow::write_parquet(all_gep_data_simulated, "all_gep_data_simulated.parquet")
}

g <- ggplot(all_gep_data_simulated) +
    geom_histogram(aes(x = DEPTH)) +
    xlim(c(0, 500)) +
    ylab("N. Events") +
    xlab("LLRG Simulated Depth") +
    ggtitle("LLRG Simulated Depth Histogram") +
    facet_wrap(. ~ Condition, scales = "free") +
    theme_bw()

ggsave("gep_simulated_hist.pdf", g, width = 15, height = 12)

g <- ggplot(all_gep_data_simulated) +
    geom_density_ridges_gradient(aes(
        x = DEPTH + 1,
        y = Condition
    )) +
    ylab("density") +
    xlim(c(0, 200)) +
    theme_ridges()

ggsave("gep_simulated_ridges.pdf", g, width = 10, height = 8)

all_gep_data_with_gene_id <- all_gep_data_simulated %>%
    dplyr::inner_join(gene_id_transcript_id, by = "TRANSCRIPT_ID") %>%
    dplyr::select(GENE_ID, TRANSCRIPT_ID, DEPTH, Condition)

if (file.exists("gep_inside_gene_isoform_level_variation_simulated.parquet")) {
    gep_inside_gene_isoform_level_variation_simulated <- arrow::read_parquet(
        "gep_inside_gene_isoform_level_variation_simulated.parquet"
    )
} else {
    gep_inside_gene_isoform_level_variation_simulated <- data.frame()
    for (gene_id in sample(unique(all_gep_data_with_gene_id$GENE_ID), 1000)) {
        for (condition in unique(all_gep_data_with_gene_id$Condition)) {
            this_gep_data_with_gene_id <- all_gep_data_with_gene_id %>%
                dplyr::filter(GENE_ID == gene_id, Condition == condition)
            this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
                dplyr::cross_join(this_gep_data_with_gene_id)
            gep_inside_gene_isoform_level_variation_simulated <- rbind(
                gep_inside_gene_isoform_level_variation_simulated,
                data.frame(
                    var = this_gep_data_with_gene_id$DEPTH.x / this_gep_data_with_gene_id$DEPTH.y,
                    Condition = condition
                )
            )
        }
    }
    arrow::write_parquet(
        gep_inside_gene_isoform_level_variation_simulated,
        "gep_inside_gene_isoform_level_variation_simulated.parquet"
    )
}

g <- ggplot(gep_inside_gene_isoform_level_variation_simulated, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Log 10 Fold Change)", limit=c(0, 5)) +
    ggtitle("Distribution of Variation Inside a Gene") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()
ggsave("gep_inside_gene_isoform_level_variation_simulated.pdf", g, width = 15, height = 12)

if (file.exists("gep_isoform_level_variation_simulated.parquet")) {
    gep_isoform_level_variation_simulated <- arrow::read_parquet("gep_isoform_level_variation_simulated.parquet")
} else {
    gep_isoform_level_variation_simulated <- data.frame()
    for (condition in unique(all_gep_data_with_gene_id$Condition)) {
        isoform_mean_depths <- all_gep_data_with_gene_id %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::select(DEPTH)
        isoform_mean_depths_min <- isoform_mean_depths %>%
            dplyr::filter(DEPTH != 0) %>%
            dplyr::arrange(DEPTH) %>%
            head(n=100)
        isoform_mean_depths_max <- isoform_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n=100)
        isoform_mean_depths <- isoform_mean_depths_max %>%
            dplyr::cross_join(isoform_mean_depths_min) %>%
            dplyr::filter(DEPTH.x >= DEPTH.y) %>%
            dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
            dplyr::select(var)
        gep_isoform_level_variation_simulated <- rbind(
            gep_isoform_level_variation_simulated,
            data.frame(
                var = isoform_mean_depths$var,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_isoform_level_variation_simulated, "gep_isoform_level_variation_simulated.parquet")
}

g <- ggplot(gep_isoform_level_variation_simulated, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Log 10 Fold Change)", limit=c(0, 5)) +
    ggtitle("Distribution of Variation among all Isoforms") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave("gep_isoform_level_variation_simulated.pdf", g, width = 15, height = 12)

if (file.exists("gep_gene_level_variation_simulated.parquet")) {
    gep_gene_level_variation_simulated <- arrow::read_parquet("gep_gene_level_variation_simulated.parquet")
} else {
    gep_gene_level_variation_simulated <- data.frame()
    for (condition in unique(all_gep_data_with_gene_id$Condition)) {
        message(sprintf("gep_gene_level_variation_simulated %s", condition))
        gene_mean_depths <- all_gep_data_with_gene_id %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::group_by(GENE_ID) %>%
            dplyr::mutate(DEPTH = mean(DEPTH)) %>%
            dplyr::distinct() %>%
            dplyr::select(DEPTH) %>%
            dplyr::filter(DEPTH != 0)
        gene_mean_depths_min <- gene_mean_depths %>%
            dplyr::arrange(DEPTH) %>%
            head(n=100)
        gene_mean_depths_max <- gene_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n=100)
        print(summary(gene_mean_depths_max))
        print(summary(gene_mean_depths_min))
        gene_mean_depths <- gene_mean_depths_max %>%
            dplyr::cross_join(gene_mean_depths_min) %>%
            dplyr::filter(DEPTH.x >= DEPTH.y) %>%
            dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
            dplyr::select(var)
        gep_gene_level_variation_simulated <- rbind(
            gep_gene_level_variation_simulated,
            data.frame(
                var = gene_mean_depths$var,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_gene_level_variation_simulated, "gep_gene_level_variation_simulated.parquet")
}

g <- ggplot(gep_gene_level_variation_simulated, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Fold Change)", limit=c(0, 5)) +
    ggtitle("Distribution of Variation among all Genes") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave("gep_gene_level_variation_simulated.pdf", g, width = 15, height = 12)
