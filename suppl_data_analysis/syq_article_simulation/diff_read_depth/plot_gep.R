library("tidyverse")
library("parallel")

args <- commandArgs(trailingOnly = TRUE)
suffix <- args[1]

cl <- parallel::makeCluster(120)


if (suffix == "simulated") {
    actual_fns <- Sys.glob("ce11_as_2_isoform_depth_*.tsv")
    conditions <- actual_fns %>%
        stringr::str_replace("ce11_as_2_isoform_depth_", "") %>%
        stringr::str_replace(".tsv", "")
} else {
    actual_fns <- Sys.glob("ce11_as_2_isoform_depth_*.fq.stats")
    conditions <- actual_fns %>%
        stringr::str_replace("ce11_as_2_isoform_depth_", "") %>%
        stringr::str_replace(".fq.stats", "")
}

all_gene_ids_with_multiple_transcripts <- readr::read_tsv(
    "ce11_as_2.gtf.gene.tsv",
    show_col_types = FALSE,
    col_types = c(
        GENE_ID = col_character(),
        TRANSCRIPT_NUMBER = col_integer(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        MAPPABLE_LENGTH = col_integer()
    )
) %>%
    dplyr::filter(TRANSCRIPT_NUMBER > 1) %>%
    dplyr::select(GENE_ID)

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


all_gep_data_filename <- sprintf("all_gep_data_%s.parquet", suffix)
if (file.exists(all_gep_data_filename)) {
    all_gep_data <- arrow::read_parquet(all_gep_data_filename)
} else {
    all_gep_data <- NULL
    for (i in seq_along(conditions)) {
        message(sprintf("%d/%d -- %s", i, length(actual_fns), actual_fns[i]))

        if (suffix == "simulated") {
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
                dplyr::mutate(Condition = conditions[i])

        } else {
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
                dplyr::rename(DEPTH = SIMULATED_DEPTH) %>%
                dplyr::right_join(all_transcript_ids, by = "TRANSCRIPT_ID") %>%
                dplyr::mutate(across(where(is.numeric), replace_na, 0)) %>%
                dplyr::mutate(
                    Condition = conditions[i]
                )
        }

        if (is.null(all_gep_data)) {
            all_gep_data <- this_gep_data
        } else {
            all_gep_data <- all_gep_data %>%
                dplyr::rows_append(this_gep_data)
        }
        rm(i, this_gep_data)
        gc()
    }
    arrow::write_parquet(all_gep_data, all_gep_data_filename)
}



if (suffix != "simulated") {
    g <- ggplot(all_gep_data) +
        geom_point(aes(
            y = DEPTH,
            x = INPUT_DEPTH
        ), size = 0.2, alpha = 0.1) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        scale_y_continuous("LLRG Simulated", trans = "log10") +
        scale_x_continuous("LLRG Input", trans = "log10") +
        ggtitle("LLRG Error Summarize") +
        facet_wrap(. ~ Condition) +
        theme_bw()

    ggsave(
        sprintf("gep_ratio_%s.png", suffix),
        g, width = 15, height = 12
    )
    means <- all_gep_data %>%
        dplyr::group_by(Condition) %>%
        dplyr::summarise(
            MEAN_INPUT_DEPTH = mean(INPUT_DEPTH),
            MEAN_SIMULATED_DEPTH = mean(DEPTH)
        )

    readr::write_tsv(means, sprintf("means_%s.tsv", suffix))
}

g <- ggplot(all_gep_data) +
    geom_histogram(aes(x = DEPTH)) +
    xlim(c(0, 500)) +
    ylab("N. Events") +
    xlab("LLRG Simulated Depth") +
    ggtitle("LLRG Simulated Depth Histogram") +
    facet_wrap(. ~ Condition, scales = "free") +
    theme_bw()
ggsave(
    sprintf("gep_hist_%s.pdf", suffix),
    g, width = 15, height = 12
)

all_gep_data_with_gene_id <- all_gep_data %>%
    dplyr::inner_join(gene_id_transcript_id, by = "TRANSCRIPT_ID") %>%
    dplyr::select(GENE_ID, TRANSCRIPT_ID, DEPTH, Condition)

all_gep_data_with_gene_id_and_multiple_transcripts <- all_gep_data_with_gene_id %>%
    dplyr::inner_join(all_gene_ids_with_multiple_transcripts, by = "GENE_ID")


gep_inside_gene_isoform_level_variation_filename <- sprintf(
    "gep_inside_gene_isoform_level_variation_%s.parquet", suffix
)
if (file.exists(gep_inside_gene_isoform_level_variation_filename)) {
    gep_inside_gene_isoform_level_variation <- arrow::read_parquet(
        gep_inside_gene_isoform_level_variation_filename
    )
} else {
    parallel::clusterExport(
        cl, 
        c("all_gep_data_with_gene_id_and_multiple_transcripts")
    )
    gep_inside_gene_isoform_level_variation_l <- parallel::parLapply(
        cl,
        unique(all_gep_data_with_gene_id_and_multiple_transcripts$GENE_ID),
        function(gene_id){
            library(dplyr)
            retdf <- data.frame()
            for (condition in unique(all_gep_data_with_gene_id_and_multiple_transcripts$Condition)) {
                this_gep_data_with_gene_id <- all_gep_data_with_gene_id_and_multiple_transcripts %>%
                    dplyr::filter(GENE_ID == gene_id, Condition == condition) %>%
                    dplyr::filter(DEPTH != 0)
                this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
                    dplyr::cross_join(this_gep_data_with_gene_id) %>%
                    dplyr::filter(DEPTH.x > DEPTH.y) %>%
                    dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
                    dplyr::select(var)
                if (length(this_gep_data_with_gene_id$var) != 0){
                    retdf <- rbind(
                        retdf,
                        data.frame(
                            var = this_gep_data_with_gene_id$var,
                            Condition = condition
                        )
                    )
                }
            }
            return(retdf)
        }
    )
    gep_inside_gene_isoform_level_variation <- Reduce(
        rbind,
        gep_inside_gene_isoform_level_variation_l
    )
    arrow::write_parquet(
        gep_inside_gene_isoform_level_variation,
        gep_inside_gene_isoform_level_variation_filename
    )
}

g <- ggplot(gep_inside_gene_isoform_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Log 10 Fold Change)", limit = c(0, 5)) +
    scale_y_continuous(limit = c(0, 600)) +
    ggtitle("Distribution of Variation Inside a Gene") +
    facet_wrap(. ~ Condition, scales="free_y") +
    theme_bw()
ggsave(
    sprintf("gep_inside_gene_isoform_level_variation_%s.pdf", suffix),
    g, width = 15, height = 12
)

gep_isoform_level_variation_filename <- sprintf(
    "gep_isoform_level_variation_%s.parquet", suffix
)
if (file.exists(gep_isoform_level_variation_filename)) {
    gep_isoform_level_variation <- arrow::read_parquet(gep_isoform_level_variation_filename)
} else {
    gep_isoform_level_variation <- data.frame()
    for (condition in unique(all_gep_data$Condition)) {
        isoform_mean_depths <- all_gep_data %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::select(DEPTH)
        isoform_mean_depths_min <- isoform_mean_depths %>%
            dplyr::filter(DEPTH != 0) %>%
            dplyr::arrange(DEPTH) %>%
            head(n = 100)
        isoform_mean_depths_max <- isoform_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n = 100)
        isoform_mean_depths <- isoform_mean_depths_max %>%
            dplyr::cross_join(isoform_mean_depths_min) %>%
            dplyr::filter(DEPTH.x >= DEPTH.y) %>%
            dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
            dplyr::select(var)
        gep_isoform_level_variation <- rbind(
            gep_isoform_level_variation,
            data.frame(
                var = isoform_mean_depths$var,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_isoform_level_variation, gep_isoform_level_variation_filename)
}

g <- ggplot(gep_isoform_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Log 10 Fold Change)", limit = c(0, 5)) +
    ggtitle("Distribution of Variation among all Isoforms") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave(
    sprintf("gep_isoform_level_variation_%s.pdf", suffix),
    g, width = 15, height = 12
)

gep_gene_level_variation_filename <- sprintf(
    "gep_gene_level_variation_%s.parquet", suffix
)
if (file.exists(gep_gene_level_variation_filename)) {
    gep_gene_level_variation <- arrow::read_parquet(gep_gene_level_variation_filename)
} else {
    gep_gene_level_variation <- data.frame()
    for (condition in unique(all_gep_data_with_gene_id$Condition)) {
        message(sprintf("gep_gene_level_variation %s", condition))
        gene_mean_depths <- all_gep_data_with_gene_id %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::group_by(GENE_ID) %>%
            dplyr::mutate(DEPTH = mean(DEPTH)) %>%
            dplyr::distinct() %>%
            dplyr::select(DEPTH) %>%
            dplyr::filter(DEPTH != 0)
        gene_mean_depths_min <- gene_mean_depths %>%
            dplyr::arrange(DEPTH) %>%
            head(n = 100)
        gene_mean_depths_max <- gene_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n = 100)
        gene_mean_depths <- gene_mean_depths_max %>%
            dplyr::cross_join(gene_mean_depths_min) %>%
            dplyr::filter(DEPTH.x >= DEPTH.y) %>%
            dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
            dplyr::select(var)
        gep_gene_level_variation <- rbind(
            gep_gene_level_variation,
            data.frame(
                var = gene_mean_depths$var,
                Condition = condition
            )
        )
    }
    arrow::write_parquet(gep_gene_level_variation, gep_gene_level_variation_filename)
}

g <- ggplot(gep_gene_level_variation, aes(x = abs(log10(var)))) +
    geom_histogram() +
    scale_x_continuous("Variance (Fold Change)", limit = c(0, 5)) +
    ggtitle("Distribution of Variation among all Genes") +
    facet_wrap(. ~ Condition, scales = "free_y") +
    theme_bw()

ggsave(
    sprintf("gep_gene_level_variation_%s.pdf", suffix),
    g, width = 15, height = 12
)
