library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")
library("DESeq2")
library("ComplexUpset")

load("dge_df_list_prepared.rd.xz")

for (i in seq_along(df_list$df)) {
    gene_covered_df <- df_list$df[[i]] %>%
        dplyr::summarise_if(is.numeric, function(x) sum(x > 0)) %>%
        t() %>%
        as_tibble(rownames = "condition", .name_repair = function(names) { "value" }) %>%
        dplyr::inner_join(
            df_list$exp,
            by = "condition"
        )
    g <- ggplot(gene_covered_df) +
        geom_boxplot(
            aes(
                x = value,
                y = sprintf("%s_%s", SIMULATOR, MODE)
            )
        ) +
        theme_bw() +
        facet_wrap(DGEID ~ DIUID) +
        ggtitle(sprintf("Number of %s covered in %s", df_list$level[[i]], df_list$stat[[i]])) +
        xlab("Sample Name") +
        ylab(sprintf("No. of %s Covered", df_list$level[[i]]))
    ggsave(
        sprintf("dge_%s_level_%s_coverage_plot.pdf", df_list$level[[i]], df_list$stat[[i]]),
        g, width = 8, height = 5
    )
}
for (i in seq_along(df_list$df)) {
    total_coverage_df <- df_list$df[[i]] %>%
        dplyr::summarise_if(is.numeric, sum) %>%
        t() %>%
        as_tibble(rownames = "condition", .name_repair = function(names) { "value" }) %>%
        dplyr::inner_join(
            df_list$exp,
            by = "condition"
        )
    g <- ggplot(total_coverage_df) +
        geom_boxplot(
            aes(
                x = value,
                y = sprintf("%s_%s", SIMULATOR, MODE)
            )
        ) +
        theme_bw() +
        facet_wrap(DGEID ~ DIUID) +
        ggtitle(sprintf("Sum of %s covered in %s", df_list$level[[i]], df_list$stat[[i]])) +
        xlab("Sample Name") +
        ylab(sprintf("Sum of %s Covered", df_list$level[[i]]))
    ggsave(
        sprintf("dge_%s_level_%s_total_coverage_plot.pdf", df_list$level[[i]], df_list$stat[[i]]),
        g, width = 8, height = 5
    )
}

gc()

dds_fc_gene_level <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared$
        all_fc_data_gene_level,
    colData = df_list$exp_prepared,
    design = ~DGEID
)
dds_fc_gene_level_vst_normalized_model <- vst(dds_fc_gene_level, blind = FALSE)
dds_fc_gene_level_vst_normalized_data <- dds_fc_gene_level_vst_normalized_model %>%
    assay() %>%
    tibble::as_tibble(rownames = "GENE_ID")

dds_fc_gene_level <- estimateSizeFactors(dds_fc_gene_level)

dds_fc_gene_level_log2p1_normalized_model <- DESeqTransform(SummarizedExperiment(
    log2(DESeq2::counts(dds_fc_gene_level, normalized = TRUE) + 1),
    colData = colData(dds_fc_gene_level)
))
dds_fc_gene_level_log2p1_normalized_data <- dds_fc_gene_level_log2p1_normalized_model %>%
    assay() %>%
    tibble::as_tibble(rownames = "GENE_ID")

dds_fc_gene_level_vst_normalized_data_mean <- dds_fc_gene_level_vst_normalized_data %>%
    dplyr::transmute(
        dge1 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_vst_normalized_data,
            tidyselect::contains("dge1")
        )),
        dge2 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_vst_normalized_data,
            tidyselect::contains("dge2")
        ))
    )
dds_fc_gene_level_log2p1_normalized_data_mean <- dds_fc_gene_level_log2p1_normalized_data %>%
    dplyr::transmute(
        dge1 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_log2p1_normalized_data,
            tidyselect::contains("dge1")
        )),
        dge2 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_log2p1_normalized_data,
            tidyselect::contains("dge2")
        ))
    )

dds_fc_gene_level_normalization_df <- dplyr::bind_rows(
    as_tibble(dds_fc_gene_level_log2p1_normalized_data_mean[, 1:2]) %>%
        dplyr::mutate(normalization_type = "log2p1"),
    as_tibble(dds_fc_gene_level_vst_normalized_data_mean[, 1:2]) %>%
        dplyr::mutate(normalization_type = "vst")
)
colnames(dds_fc_gene_level_normalization_df)[1:2] <- c("dge1", "dge2")
p <- ggplot(
    dds_fc_gene_level_normalization_df,
    aes(x = dge1, y = dge2)
) +
    geom_hex(bins = 80) +
    coord_fixed() +
    facet_grid(. ~ normalization_type) +
    labs(xlab = "X",
         ylab = "Y",
         title = "Normalization method"
    )
ggsave(
    sprintf("dge_fc_gene_level_normalized.pdf"),
    p, width = 8, height = 5
)

p <- DESeq2::plotPCA(
    dds_fc_gene_level_vst_normalized_model,
    intgroup = "DGEID"
)
ggsave(
    sprintf("dge_fc_gene_level_vst_pca.pdf"),
    p, width = 5, height = 5
)

dds_fc_gene_level <- DESeq(dds_fc_gene_level)

dds_fc_gene_level_normalized_counts <- dds_fc_gene_level %>%
    DESeq2::counts(, normalized = TRUE) %>%
    as_tibble(rownames = "GENE_ID")
dds_fc_gene_level_basemean <- dds_fc_gene_level_normalized_counts %>%
    dplyr::transmute(
        dge1 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_normalized_counts,
            tidyselect::contains("dge1")
        )),
        dge2 = base::rowMeans(dplyr::select(
            dds_fc_gene_level_normalized_counts,
            tidyselect::contains("dge2")
        ))
    )
arrow::write_parquet(dds_fc_gene_level_basemean, "dge_fc_gene_level_base_mean.parquet")
dds_fc_gene_level_dge_result <- dds_fc_gene_level %>%
    DESeq2::results(lfcThreshold = 2, alpha = 0.05) %>%
    tibble::as_tibble(, rownames = "GENE_ID")
arrow::write_parquet(dds_fc_gene_level_dge_result, "dge_fc_gene_level_result.parquet")

pdf(sprintf("dge_fc_gene_level_ma.pdf"), width = 8, height = 5)
plotMA(dds_fc_gene_level)
dev.off()

dds_fc_gene_level_dge_result_filtered <- dds_fc_gene_level_dge_result %>%
    dplyr::filter(
        abs(log2FoldChange) >= 2,
        padj < 0.05
    )
arrow::write_parquet(dds_fc_gene_level_dge_result_filtered, "dge_fc_gene_level_result_filtered.parquet")

dds_fc_gene_level_pheatmap_data <- dds_fc_gene_level_vst_normalized_data %>%
    as.data.frame()
row.names(dds_fc_gene_level_pheatmap_data) <- dds_fc_gene_level_vst_normalized_data$GENE_ID
dds_fc_gene_level_pheatmap_data <- dds_fc_gene_level_pheatmap_data %>%
    dplyr::select(!(GENE_ID)) %>%
    dplyr::mutate(var = matrixStats::rowVars(as.matrix(.))) %>%
    dplyr::arrange(desc(var)) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::select(!var)

repless_conditions <- (df_list$exp_prepared %>%
    dplyr::select(!REPID) %>%
    dplyr::mutate(repless_condition = paste(SIMULATOR, SEQUENCER, MODE, DGEID, DIUID, sep = "_")))$repless_condition

dds_fc_gene_level_pheatmap_data_means <- data.frame(
    row.names = row.names(dds_fc_gene_level_pheatmap_data)
)
for (rep_name in repless_conditions) {
    dds_fc_gene_level_pheatmap_data_means[[rep_name]] <- dds_fc_gene_level_pheatmap_data %>%
        dplyr::select(tidyselect::contains(rep_name)) %>%
        base::rowMeans()
    rm(rep_name)
}

pheatmap(
    dds_fc_gene_level_pheatmap_data_means,
    height = 10,
    width = 8,
    filename = "dge_fc_gene_level_heatmap.pdf"
)

#' Detect DGE Overlapping Rate
dds_fc_gene_level_dge_result <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared$
        all_fc_data_gene_level,
    colData = df_list$exp_prepared,
    design = ~DGEID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "GENE_ID") %>%
    dplyr::transmute(
        GENE_ID = GENE_ID,
        is_dge_fc = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_fq_gene_level_dge_result <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared$
        all_fq_stats_gene_level,
    colData = df_list$exp_prepared,
    design = ~DGEID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "GENE_ID") %>%
    dplyr::transmute(
        GENE_ID = GENE_ID,
        is_dge_fq = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_gene_level_dge_result <- tibble::tibble(
    GENE_ID = Reduce(
        base::union,
        c(
            dds_fc_gene_level_dge_result$GENE_ID,
            dds_fq_gene_level_dge_result$GENE_ID,
            as_transcripts_data$GENE_ID,
            reference_transcripts_data$GENE_ID
        )
    )
) %>%
    dplyr::left_join(dds_fc_gene_level_dge_result, by = "GENE_ID") %>%
    dplyr::left_join(dds_fq_gene_level_dge_result, by = "GENE_ID")

pdf("dge_gene_level_upset.pdf", width = 8, height = 5)
upset(
    dds_gene_level_dge_result,
    intersect = c("is_dge_fq", "is_dge_fc")
)
dev.off()


dds_fc_isoform_level_dge_result <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared$
        all_fc_data_isoform_level,
    colData = df_list$exp_prepared,
    design = ~DGEID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "TRANSCRIPT_ID") %>%
    dplyr::transmute(
        TRANSCRIPT_ID = TRANSCRIPT_ID,
        is_dge_fc = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_fq_isoform_level_dge_result <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared$
        all_fq_stats_isoform_level,
    colData = df_list$exp_prepared,
    design = ~DGEID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "TRANSCRIPT_ID") %>%
    dplyr::transmute(
        TRANSCRIPT_ID = TRANSCRIPT_ID,
        is_dge_fq = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_isoform_level_dge_result <- tibble::tibble(
    TRANSCRIPT_ID = Reduce(
        base::union,
        c(
            dds_fc_isoform_level_dge_result$TRANSCRIPT_ID,
            dds_fq_isoform_level_dge_result$GENE_ID,
            as_transcripts_data$TRANSCRIPT_ID,
            reference_transcripts_data$TRANSCRIPT_ID
        )
    )
) %>%
    dplyr::left_join(dds_fc_isoform_level_dge_result, by = "TRANSCRIPT_ID") %>%
    dplyr::left_join(dds_fq_isoform_level_dge_result, by = "TRANSCRIPT_ID")

pdf("dge_isoform_level_upset.pdf", width = 8, height = 5)
upset(
    dds_isoform_level_dge_result,
    intersect = c("is_dge_fq", "is_dge_fc")
)
dev.off()

rm(list = grep("dds", ls(), value = TRUE))
gc()

#' DGE accross DIU analysis


dds_fc_gene_level_dge_result_dge1_only <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared_dge1_only$
        all_fc_data_gene_level,
    colData = df_list$exp_prepared_dge1_only,
    design = ~DIUID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "GENE_ID") %>%
    dplyr::transmute(
        GENE_ID = GENE_ID,
        is_dge_fc = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_fq_gene_level_dge_result_dge1_only <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared_dge1_only$
        all_fq_stats_gene_level,
    colData = df_list$exp_prepared_dge1_only,
    design = ~DIUID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "GENE_ID") %>%
    dplyr::transmute(
        GENE_ID = GENE_ID,
        is_dge_fq = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )


as_transcripts_data <- readr::read_tsv(
    "ce11.ncbiRefSeq_as.chr1.gtf.transcripts.tsv",
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
)
reference_transcripts_data <- readr::read_tsv(
    "ce11.ncbiRefSeq.chr1.gtf.transcripts.tsv",
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
)
dds_gene_level_dge_result_dge1_only <- tibble::tibble(
    GENE_ID = Reduce(
        base::union,
        c(
            dds_fc_gene_level_dge_result_dge1_only$GENE_ID,
            dds_fq_gene_level_dge_result_dge1_only$GENE_ID,
            as_transcripts_data$GENE_ID,
            reference_transcripts_data$GENE_ID
        )
    )
) %>%
    dplyr::left_join(dds_fc_gene_level_dge_result_dge1_only, by = "GENE_ID") %>%
    dplyr::left_join(dds_fq_gene_level_dge_result_dge1_only, by = "GENE_ID")

pdf("dge_gene_level_upset_dge1_only.pdf", width = 8, height = 5)
upset(
    dds_gene_level_dge_result_dge1_only,
    intersect = c("is_dge_fq", "is_dge_fc")
)
dev.off()


dds_fc_isoform_level_dge_result_dge1_only <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared_dge1_only$
        all_fc_data_isoform_level,
    colData = df_list$exp_prepared_dge1_only,
    design = ~DIUID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "TRANSCRIPT_ID") %>%
    dplyr::transmute(
        TRANSCRIPT_ID = TRANSCRIPT_ID,
        is_dge_fc = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_fq_isoform_level_dge_result_dge1_only <- DESeqDataSetFromMatrix(
    countData = df_list$
        df_filtered_prepared_dge1_only$
        all_fq_stats_isoform_level,
    colData = df_list$exp_prepared_dge1_only,
    design = ~DIUID
) %>%
    DESeq() %>%
    DESeq2::results() %>%
    tibble::as_tibble(rownames = "TRANSCRIPT_ID") %>%
    dplyr::transmute(
        TRANSCRIPT_ID = TRANSCRIPT_ID,
        is_dge_fq = (abs(log2FoldChange) >= 2 & padj < 0.05)
    )

dds_isoform_level_dge_result_dge1_only <- tibble::tibble(
    TRANSCRIPT_ID = Reduce(
        base::union,
        c(
            dds_fc_isoform_level_dge_result_dge1_only$TRANSCRIPT_ID,
            dds_fq_isoform_level_dge_result_dge1_only$GENE_ID,
            as_transcripts_data$TRANSCRIPT_ID,
            reference_transcripts_data$TRANSCRIPT_ID
        )
    )
) %>%
    dplyr::left_join(dds_fc_isoform_level_dge_result_dge1_only, by = "TRANSCRIPT_ID") %>%
    dplyr::left_join(dds_fq_isoform_level_dge_result_dge1_only, by = "TRANSCRIPT_ID")

pdf("dge_isoform_level_upset_dge1_only.pdf", width = 8, height = 5)
upset(
    dds_isoform_level_dge_result_dge1_only,
    intersect = c("is_dge_fq", "is_dge_fc")
)
dev.off()


rm(list = grep("dds", ls(), value = TRUE))
gc()
