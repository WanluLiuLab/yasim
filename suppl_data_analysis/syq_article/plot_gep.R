library("tidyverse")
library("parallel")

cl <- parallel::makeCluster(40)


actual_fns <- Sys.glob("real_expression/*.parquet")
conditions <- actual_fns %>%
    stringr::str_replace("real_expression/", "") %>%
    stringr::str_replace(".fastq.trans.depth.filtered.parquet", "") %>%
    stringr::str_replace(".fastq.gz.trans.depth.filtered.parquet", "")

all_gep_data_filename <- sprintf("all_gep_data.parquet")
if (file.exists(all_gep_data_filename)) {
    all_gep_data <- arrow::read_parquet(all_gep_data_filename)
} else {
    all_gep_data <- NULL
    for (i in seq_along(conditions)) {
        message(sprintf("%d/%d -- %s", i, length(actual_fns), actual_fns[i]))
        this_gep_data <- arrow::read_parquet(
            actual_fns[i]
        ) %>%
            dplyr::select(TRANSCRIPT_ID, GENE_ID, DEPTH) %>%
            dplyr::mutate(Condition = conditions[i])

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

all_gep_data <- all_gep_data %>%
    dplyr::filter(DEPTH >= 1)
g <- all_gep_data %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(DEPTH = mean(DEPTH)) %>%
    ggplot() +
    geom_bar(aes(x = DEPTH, y = Condition, fill = Condition), stat = "identity") +
    xlab("Depth") +
    ggtitle("Depth Bar plot") +
    theme_bw()
ggsave(
    sprintf("gep_bar.pdf"),
    g,
    width = 15, height = 12
)
quit()

# g <- ggplot(all_gep_data) +
#     geom_histogram(aes(x = DEPTH)) +
#     scale_y_continuous("N. Events", trans = "log10") +
#     xlab("Depth") +
#     ggtitle("Depth Histogram") +
#     facet_wrap(. ~ Condition, scales = "free") +
#     theme_bw()
# ggsave(
#     sprintf("gep_hist.pdf"),
#     g,
#     width = 15, height = 12
# )

gep_inside_gene_isoform_level_variation_filename <- "gep_inside_gene_isoform_level_variation.parquet"

if (file.exists(gep_inside_gene_isoform_level_variation_filename)) {
    gep_inside_gene_isoform_level_variation <- arrow::read_parquet(
        gep_inside_gene_isoform_level_variation_filename
    )
} else {
    parallel::clusterExport(
        cl,
        c("all_gep_data")
    )
    gep_inside_gene_isoform_level_variation_l <- parallel::parLapply(
        cl,
        conditions,
        function(condition) {
            library(dplyr)
            retdf <- data.frame()

            gene_ids <- (
                all_gep_data %>%
                    dplyr::filter(Condition == condition) %>%
                    dplyr::group_by(GENE_ID) %>%
                    dplyr::summarise(n = n()) %>%
                    dplyr::ungroup() %>%
                    dplyr::filter(n >= 2) %>%
                    dplyr::select(GENE_ID)
            )$GENE_ID
            # gene_ids <- sample(
            #     gene_ids,
            #     as.integer(length(gene_ids) * 0.05)
            # )
            this_gep_data <- all_gep_data %>%
                    dplyr::filter(Condition == condition) %>%
                    dplyr::filter(GENE_ID %in% gene_ids)
            for (gene_id in gene_ids) {
                this_gep_data_with_gene_id <- this_gep_data %>%
                    dplyr::filter(GENE_ID == gene_id)
                if (nrow(this_gep_data_with_gene_id) >= 2) {
                    this_gep_data_with_gene_id <- this_gep_data_with_gene_id %>%
                        dplyr::cross_join(this_gep_data_with_gene_id) %>%
                        dplyr::filter(DEPTH.x > DEPTH.y) %>%
                        dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
                        dplyr::select(var)
                    if (nrow(this_gep_data_with_gene_id) != 0) {
                        retdf <- rbind(
                            retdf,
                            data.frame(
                                var = this_gep_data_with_gene_id$var,
                                Condition = condition
                            )
                        )
                    }
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
# g <- ggplot(gep_inside_gene_isoform_level_variation, aes(x = abs(log10(var)))) +
#     geom_histogram() +
#     scale_x_continuous("Variance (Log 10 Fold Change)", limit = c(0, 5)) +
#     scale_y_continuous(limit = c(0, 600)) +
#     ggtitle("Distribution of Variation Inside a Gene") +
#     facet_wrap(. ~ Condition, scales = "free_y") +
#     theme_bw()
# ggsave(
#     "gep_inside_gene_isoform_level_variation.pdf",
#     g, width = 15, height = 12
# )
rm(gep_inside_gene_isoform_level_variation)
gc()

gep_isoform_level_variation_filename <- "gep_isoform_level_variation.parquet"
if (file.exists(gep_isoform_level_variation_filename)) {
    gep_isoform_level_variation <- arrow::read_parquet(gep_isoform_level_variation_filename)
} else {
    gep_isoform_level_variation <- data.frame()
    for (condition in unique(all_gep_data$Condition)) {
        message(sprintf("gep_isoform_level_variation %s", condition))
        isoform_mean_depths <- all_gep_data %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::select(DEPTH)
        gc()
        isoform_mean_depths_min <- isoform_mean_depths %>%
            dplyr::arrange(DEPTH) %>%
            head(n = 0.025 * nrow(isoform_mean_depths))
        isoform_mean_depths_max <- isoform_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n = 0.025 * nrow(isoform_mean_depths))
        rm(isoform_mean_depths)
        gc()
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
        gc()
    }
    arrow::write_parquet(gep_isoform_level_variation, gep_isoform_level_variation_filename)
}

# g <- ggplot(gep_isoform_level_variation, aes(x = abs(log10(var)))) +
#     geom_histogram() +
#     scale_x_continuous("Variance (Log 10 Fold Change)", limit = c(0, 5)) +
#     ggtitle("Distribution of Variation among all Isoforms") +
#     facet_wrap(. ~ Condition, scales = "free_y") +
#     theme_bw()

# ggsave(
#     "gep_isoform_level_variation.pdf",
#     g, width = 15, height = 12
# )
rm(gep_isoform_level_variation)
gc()

gep_gene_level_variation_filename <- "gep_gene_level_variation.parquet"

if (file.exists(gep_gene_level_variation_filename)) {
    gep_gene_level_variation <- arrow::read_parquet(gep_gene_level_variation_filename)
} else {
    gep_gene_level_variation <- data.frame()
    for (condition in unique(all_gep_data$Condition)) {
        message(sprintf("gep_gene_level_variation %s", condition))
        gene_mean_depths <- all_gep_data %>%
            dplyr::filter(Condition == condition) %>%
            dplyr::group_by(GENE_ID) %>%
            dplyr::mutate(DEPTH = mean(DEPTH)) %>%
            dplyr::ungroup() %>%
            dplyr::distinct() %>%
            dplyr::select(DEPTH)
        gc()
        gene_mean_depths_min <- gene_mean_depths %>%
            dplyr::arrange(DEPTH) %>%
            head(n = 0.025 * nrow(gene_mean_depths))
        gene_mean_depths_max <- gene_mean_depths %>%
            dplyr::arrange(desc(DEPTH)) %>%
            head(n = 0.025 * nrow(gene_mean_depths))
        rm(gene_mean_depths)
        gc()
        gene_mean_depths <- gene_mean_depths_max %>%
            dplyr::cross_join(gene_mean_depths_min) %>%
            dplyr::filter(DEPTH.x >= DEPTH.y) %>%
            dplyr::mutate(var = DEPTH.x / DEPTH.y) %>%
            dplyr::select(var)
        gc()
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

# g <- ggplot(gep_gene_level_variation, aes(x = abs(log10(var)))) +
#     geom_histogram() +
#     scale_x_continuous("Variance (Fold Change)", limit = c(0, 5)) +
#     ggtitle("Distribution of Variation among all Genes") +
#     facet_wrap(. ~ Condition, scales = "free_y") +
#     theme_bw()

# ggsave(
#     "gep_gene_level_variation.pdf",
#     g,
#     width = 15, height = 12
# )
