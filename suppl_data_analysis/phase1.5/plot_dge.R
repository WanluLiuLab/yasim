library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")
library("DESeq2")
library("ComplexUpset")


#' Merging Data

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

all_fc_data_gene_level <- NULL
all_fc_data_isoform_level <- NULL
experiment_design <- NULL
fns <- Sys.glob("ce11_*.fq.gz.bam.fc.tsv")
conditions <- fns %>%
    stringr::str_replace(".fq.gz.bam.fc.tsv", "") %>%
    stringr::str_replace("ce11_", "")

for (i in seq_along(fns)) {
    fc_data_fn <- fns[i]
    condition <- conditions[i]
    this_fc_data <- readr::read_tsv(
        fc_data_fn,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            Chr = col_character(),
            Start = col_character(),
            End = col_character(),
            Strand = col_character(),
            Length = col_integer(),
            NumReads = col_integer(),
        ),
        col_names = c(
            "TRANSCRIPT_ID",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
            "NumReads"
        ),
        comment = "#"
    ) %>%
        dplyr::select(TRANSCRIPT_ID, NumReads) %>%
        tidyr::drop_na()
    this_fc_data_isoform_level <- this_fc_data %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(condition) := NumReads
        )
    this_fc_data_gene_level <- this_fc_data %>%
        dplyr::inner_join(
            reference_transcripts_data,
            by = "TRANSCRIPT_ID"
        ) %>%
        dplyr::group_by(GENE_ID) %>%
        dplyr::summarise(NumReads = sum(NumReads)) %>%
        dplyr::transmute(
            GENE_ID = GENE_ID,
            !!rlang::sym(condition) := NumReads
        )
    if (is.null(all_fc_data_isoform_level)) {
        all_fc_data_isoform_level <- this_fc_data_isoform_level %>%
            dplyr::select(TRANSCRIPT_ID)
    }
    if (is.null(all_fc_data_gene_level)) {
        all_fc_data_gene_level <- this_fc_data_gene_level %>%
            dplyr::select(GENE_ID)
    }
    all_fc_data_isoform_level <- all_fc_data_isoform_level %>%
        dplyr::inner_join(this_fc_data_isoform_level, by = "TRANSCRIPT_ID")
    all_fc_data_gene_level <- all_fc_data_gene_level %>%
        dplyr::inner_join(this_fc_data_gene_level, by = "GENE_ID")

    condition_break <- strsplit(condition, "_")
    this_experiment_design <- data.frame(
        SIMULATOR = condition_break[[1]][1],
        MODE = condition_break[[1]][2],
        DGEID = condition_break[[1]][3],
        DIUID = condition_break[[1]][4],
        REPID = condition_break[[1]][5],
        condition = condition
    )
    if (is.null(experiment_design)) {
        experiment_design <- this_experiment_design
    } else {
        experiment_design <- experiment_design %>%
            dplyr::rows_append(
                this_experiment_design
            )
    }
    message(sprintf("Processing %d/%d", i, length(fns)))
    rm(
        condition,
        condition_break,
        this_fc_data,
        this_fc_data_isoform_level,
        this_fc_data_gene_level,
        fc_data_fn,
        this_experiment_design,
        i
    )
    gc()
}

arrow::write_parquet(all_fc_data_gene_level, "dge_fc_data_gene_level.parquet")
arrow::write_parquet(all_fc_data_isoform_level, "dge_fc_data_isoform_level.parquet")
arrow::write_parquet(experiment_design, "dge_experiment_design.parquet")

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

all_fq_stats_gene_level <- NULL
all_fq_stats_isoform_level <- NULL
fns <- Sys.glob("ce11_*.fq.stats")
conditions <- fns %>%
    stringr::str_replace(".fq.stats", "") %>%
    stringr::str_replace("ce11_", "")

for (i in seq_along(fns)) {
    condition <- conditions[i]
    this_fq_stats_data <- readr::read_tsv(
        fns[i],
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            INPUT_DEPTH = col_integer(),
            SIMULATED_N_OF_READS = col_integer()
        ),
        comment = "#"
    )
    this_fq_stats_data_isoform_level <- this_fq_stats_data %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(condition) := SIMULATED_N_OF_READS
        )
    this_fq_stats_data_gene_level <- this_fq_stats_data %>%
        dplyr::inner_join(
            reference_transcripts_data,
            by = "TRANSCRIPT_ID"
        ) %>%
        dplyr::group_by(GENE_ID) %>%
        dplyr::summarise(SIMULATED_N_OF_READS = sum(SIMULATED_N_OF_READS)) %>%
        dplyr::transmute(
            GENE_ID = GENE_ID,
            !!rlang::sym(condition) := SIMULATED_N_OF_READS
        )

    if (is.null(all_fq_stats_isoform_level)) {
        all_fq_stats_isoform_level <- this_fq_stats_data_isoform_level
    } else {
        all_fq_stats_isoform_level <- all_fq_stats_isoform_level %>%
            dplyr::inner_join(
                this_fq_stats_data_isoform_level,
                by = "TRANSCRIPT_ID"
            )
    }
    if (is.null(all_fq_stats_gene_level)) {
        all_fq_stats_gene_level <- this_fq_stats_data_gene_level
    } else {
        all_fq_stats_gene_level <- all_fq_stats_gene_level %>%
            dplyr::inner_join(
                this_fq_stats_data_gene_level,
                by = "GENE_ID"
            )
    }
    message(sprintf("Processing %d/%d -- %s", i, length(fns), fns[i]))
    rm(
        this_fq_stats_data,
        this_fq_stats_data_isoform_level,
        this_fq_stats_data_gene_level,
        condition,
        i
    )
    gc()
}

arrow::write_parquet(all_fq_stats_gene_level, "dge_fq_stats_gene_level.parquet")
arrow::write_parquet(all_fq_stats_isoform_level, "dge_fq_stats_isoform_level.parquet")

#' Basic Stats

df_list <- list(
    df = list(
        all_fq_stats_gene_level = all_fq_stats_gene_level,
        all_fq_stats_isoform_level = all_fq_stats_isoform_level,
        all_fc_data_gene_level = all_fc_data_gene_level,
        all_fc_data_isoform_level = all_fc_data_isoform_level
    ),
    level = c("gene", "isoform", "gene", "isoform"),
    stat = c("YASIM", "YASIM", "FC", "FC"),
    exp = experiment_design
)

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
            experiment_design,
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

#' Call DGE
filter_genes_isoforms <- function(df) {
    df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(total = sum(c_across(tidyselect::where(is.numeric)))) %>%
        dplyr::filter(total > 2 * length(df)) %>%
        dplyr::select(!(total))
}

df_list[["df_filtered"]] <- list(
    all_fq_stats_gene_level = filter_genes_isoforms(df_list$df$all_fq_stats_gene_level),
    all_fq_stats_isoform_level = filter_genes_isoforms(df_list$df$all_fq_stats_isoform_level),
    all_fc_data_gene_level = filter_genes_isoforms(df_list$df$all_fc_data_gene_level),
    all_fc_data_isoform_level = filter_genes_isoforms(df_list$df$all_fc_data_isoform_level)
)

prepare_tibble_for_dge <- function(df, row_name) {
    df <- as.data.frame(df)
    rownames(df) <- df[[row_name]]
    df[[row_name]] <- NULL
    return(df)
}

df_list[["df_filtered_prepared"]] <- list(
    all_fq_stats_gene_level = prepare_tibble_for_dge(
        df_list$df_filtered$all_fq_stats_gene_level, "GENE_ID"
    ),
    all_fq_stats_isoform_level = prepare_tibble_for_dge(
        df_list$df_filtered$all_fq_stats_isoform_level, "TRANSCRIPT_ID"
    ),
    all_fc_data_gene_level = prepare_tibble_for_dge(
        df_list$df_filtered$all_fc_data_gene_level, "GENE_ID"
    ),
    all_fc_data_isoform_level = prepare_tibble_for_dge(
        df_list$df_filtered$all_fc_data_isoform_level, "TRANSCRIPT_ID"
    )
)
df_list[["exp_prepared"]] <- prepare_tibble_for_dge(
    df_list$exp, "condition"
)

rm(
    all_fc_data_gene_level,
    all_fc_data_isoform_level,
    all_fq_stats_isoform_level,
    all_fq_stats_gene_level,
    experiment_design,
    conditions,
    fns
)
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

repless_conditions <- (df_list$exp %>%
    dplyr::select(!REPID) %>%
    dplyr::mutate(repless_condition = paste(SIMULATOR, MODE, DGEID, DIUID, sep = "_")))$repless_condition

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
rm(list = grep("dds", ls(), value = TRUE))
rm(g)
gc()

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

prepare_tibble_for_dge_dge1_only <- function(df, row_name) {
    df <- df %>%
        dplyr::select(tidyr::contains("dge1")) %>%
        as.data.frame(df)
    rownames(df) <- df[[row_name]]
    df[[row_name]] <- NULL
    return(df)
}

df_list[["df_filtered_prepared_dge1_only"]] <- list(
    all_fq_stats_gene_level = filter_genes_isoforms(prepare_tibble_for_dge_dge1_only(
        df_list$df$all_fq_stats_gene_level, "GENE_ID"
    )),
    all_fq_stats_isoform_level = filter_genes_isoforms(prepare_tibble_for_dge_dge1_only(
        df_list$df$all_fq_stats_isoform_level, "TRANSCRIPT_ID"
    )),
    all_fc_data_gene_level = filter_genes_isoforms(prepare_tibble_for_dge_dge1_only(
        df_list$df$all_fc_data_gene_level, "GENE_ID"
    )),
    all_fc_data_isoform_level = filter_genes_isoforms(prepare_tibble_for_dge_dge1_only(
        df_list$df$all_fc_data_isoform_level, "TRANSCRIPT_ID"
    ))
)
df_list[["exp_prepared_dge1_only"]] <- prepare_tibble_for_dge(
    df_list$exp %>% dplyr::filter(DGEID == "dge1"), "conditions"
)

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
save(
    df_list = df_list,
    file = "dge_df_list_prepared.rd.xz",
    compress = "xz",
    compression_level = 9
)

rm(list = grep("dds", ls(), value = TRUE))
gc()
