library(tidyverse)
library(DESeq2)

full_df <- arrow::read_parquet("dge_fq_stats_isoform_level.parquet")
exp_design <- arrow::read_parquet("dge_experiment_design.parquet") %>%
    dplyr::filter(DGEID == "dge1")
reference_transcripts_data <- readr::read_tsv(
    "ce11.ncbiRefSeq_as.chr1.gtf.transcripts.tsv",
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
)

full_df <- data.frame(
    TRANSCRIPT_ID = c("G1.1", "G1.2", "G2.1", "G2.2"),
    diu1_1 = c(1, 5, 4, 4),
    diu1_2 = c(2, 4, 4, 4),
    diu2_1 = c(4, 1, 16, 16),
    diu2_2 = c(5, 2, 16, 16)
)
exp_design <- data.frame(
    condition = c("diu1_1", "diu1_2", "diu2_1", "diu2_2"),
    DIUID = c("diu1", "diu1", "diu2", "diu2")
)
reference_transcripts_data <- data.frame(
    TRANSCRIPT_ID = c("G1.1", "G1.2", "G2.1", "G2.2"),
    GENE_ID = c("G1", "G1", "G2", "G2")
)

# cnames <- grep("dge1", full_df %>% colnames(), value = TRUE)
cnames <- full_df %>% colnames()
cnames1 <- grep("diu1", cnames, value = TRUE)
cnames2 <- grep("diu2", cnames, value = TRUE)

full_df <- full_df %>%
    # dplyr::select(tidyselect::contains("dge1"), TRANSCRIPT_ID) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = prod(c_across(tidyselect::where(is.numeric)))) %>%
    dplyr::filter(total != 0) %>%
    dplyr::select(!(total))


gene_transcript_df_detected <- full_df %>%
    dplyr::select(TRANSCRIPT_ID) %>%
    dplyr::inner_join(
        reference_transcripts_data,
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::select(TRANSCRIPT_ID, GENE_ID)

gene_transcript_table <- list()
for (gene_id in unique(gene_transcript_df_detected$GENE_ID)) {
    this_transcript_ids <- gene_transcript_df_detected %>%
        dplyr::filter(GENE_ID == gene_id) %>%
        dplyr::select(TRANSCRIPT_ID) %>%
        unlist() %>%
        as.character()
    if (length(this_transcript_ids) >= 2) {
        gene_transcript_table[[gene_id]] <- this_transcript_ids
    }
    rm(gene_id, this_transcript_ids)
}

rm(gene_transcript_df_detected)

full_df_mutated <- as.data.frame(full_df)
rownames(full_df_mutated) <- full_df_mutated$TRANSCRIPT_ID
full_df_mutated$TRANSCRIPT_ID <- NULL

c_df <- data.frame()
gene_ids <- names(gene_transcript_table)


for (i in seq_along(gene_ids)) {
    gene_id <- gene_ids[i]
    print(sprintf("%s -- %d/%d", gene_id, i, length(gene_ids)))
    sel_df <- full_df_mutated[gene_transcript_table[[gene_id]],]
    sel_dds <- suppressWarnings(
        DESeq2::DESeqDataSetFromMatrix(
            countData = sel_df,
            colData = exp_design,
            design = ~DIUID
        ) %>%
            DESeq2::estimateSizeFactors()
    )
    sel_dds_nc <- DESeq2::counts(sel_dds, normalize = TRUE) %>%
        as.data.frame()
    sel_dds_nc_1 <- sel_dds_nc %>%
        dplyr::select(tidyselect::contains("diu1"))
    sel_dds_nc_2 <- sel_dds_nc %>%
        dplyr::select(tidyselect::contains("diu2"))
    for (transcript_id in row.names(sel_dds_nc)) {
        l2fc <- log2(
            mean(as.double(sel_dds_nc_1[transcript_id,])) /
                mean(as.double(sel_dds_nc_2[transcript_id,]))
        )
        basemean <- mean(as.double(sel_dds_nc[transcript_id,]))
        pv <- suppressWarnings(wilcox.test(
            as.double(sel_dds_nc_1[transcript_id,]),
            as.double(sel_dds_nc_2[transcript_id,])
        )$p.value)
        c_df <- rbind(
            c_df,
            data.frame(
                GENE_ID = gene_id,
                TRANSCRIPT_ID = transcript_id,
                pv = pv,
                basemean = basemean,
                l2fc = l2fc
            )
        )
    }
}

ggplot(c_df) +
    geom_point(aes(
        x = l2fc,
        y = basemean,
        color = pv < 0.05 & abs(l2fc) > 2
    )) +
    geom_label(aes(
        x = l2fc,
        y = basemean,
        label = TRANSCRIPT_ID
    )) +
    scale_color_manual(
        values = c("black", "red")
    ) +
    theme_bw()

ggplot(c_df) +
    geom_point(aes(
        x = l2fc,
        y = log10(pv),
        color = pv < 0.05 & abs(l2fc) > 2
    )) +
    scale_color_manual(
        values = c("black", "red")
    ) +
    theme_bw()
