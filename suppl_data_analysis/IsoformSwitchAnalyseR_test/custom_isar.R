library(tidyverse)
library(DESeq2)

full_df <- arrow::read_parquet("dge_fq_stats_isoform_level.parquet")
exp_design <- arrow::read_parquet("dge_experiment_design.parquet") %>%
    dplyr::filter(DGEID=="dge1")

cnames <- grep("dge1", full_df %>% colnames(), value = TRUE)
cnames1 <- grep("diu1", cnames, value = TRUE)
cnames2 <- grep("diu2", cnames, value = TRUE)

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

gene_transcript_df_detected <- full_df %>%
    dplyr::select(TRANSCRIPT_ID) %>%
    dplyr::inner_join(
        reference_transcripts_data,
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::select(TRANSCRIPT_ID, GENE_ID)

full_df <- full_df %>%
    dplyr::inner_join(
        gene_transcript_df_detected,
        by = "TRANSCRIPT_ID"
    ) %>%
    dplyr::select(tidyselect::contains("dge1"), GENE_ID, TRANSCRIPT_ID) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(total = prod(c_across(tidyselect::where(is.numeric)))) %>%
        dplyr::filter(total != 0) %>%
        dplyr::select(!(total))

gene_ids <- unique(full_df$GENE_ID)

c_df <- data.frame()

i <- 0
for (gene_id in gene_ids) {
    print(sprintf("%s -- %d/%d", gene_id, i, length(gene_ids)))
    i <- i + 1
    sel_df <- full_df %>%
        dplyr::filter(GENE_ID == gene_id) %>%
        as.data.frame()
    row.names(sel_df) <- sel_df$TRANSCRIPT_ID
    sel_df$TRANSCRIPT_ID <- NULL
    sel_df$GENE_ID <- NULL
    # print(sel_df)
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
    for (transcript_id in row.names(sel_dds_nc)){
        l2fc <- log2(
            mean(as.double(sel_dds_nc_1[transcript_id, ])) /
                mean(as.double(sel_dds_nc_2[transcript_id, ]))
        )
        basemean <- mean(as.double(sel_dds_nc[transcript_id, ]))
        pv <- suppressWarnings(wilcox.test(
            as.double(sel_dds_nc_1[transcript_id, ]),
            as.double(sel_dds_nc_2[transcript_id, ])
        )$p.value)
        c_df <- rbind(
            c_df,
            data.frame(
                GENE_ID = gene_id,
                TRANSCRIPT_ID = transcript_id,
                pv = pv,
                basemean=basemean,
                l2fc = l2fc
            )
        )
    }
}

ggplot(c_df) +
    geom_point(aes(
        x=l2fc,
        y=basemean,
        color=pv < 0.05 & abs(l2fc) > 2
    ))
