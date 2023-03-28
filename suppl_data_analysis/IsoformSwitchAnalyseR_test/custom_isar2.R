library(tidyverse)
library(DESeq2)

library(parallel)

cl <- parallel::makeCluster(parallel::detectCores())

fine2a <- readr::read_tsv(
    "FINE2a.txt",
    comment = "#"
) %>%
    dplyr::transmute(
        TRANSCRIPT_ID = Geneid,
        FINE2a=FINE2a.bam
    )

fine2b <- readr::read_tsv(
    "FINE2b.txt",
    comment = "#"
) %>%
    dplyr::transmute(
        TRANSCRIPT_ID = Geneid,
        FINE2b=FINE2b.bam
    )

tesr7a <- readr::read_tsv(
    "TesR7A.txt",
    comment = "#"
) %>%
    dplyr::transmute(
        TRANSCRIPT_ID = Geneid,
        TesR7A=TesR7A.bam
    )

tesr7b <- readr::read_tsv(
    "TesR7B.txt",
    comment = "#"
) %>%
    dplyr::transmute(
        TRANSCRIPT_ID = Geneid,
        TesR7B=TesR7B.bam
    )

full_df <- fine2a %>%
    dplyr::inner_join(
            fine2b,
        by="TRANSCRIPT_ID"
    ) %>%
        dplyr::inner_join(
            tesr7a,
        by="TRANSCRIPT_ID"
    ) %>%
        dplyr::inner_join(
            tesr7b,
        by="TRANSCRIPT_ID"
    )
arrow::write_parquet(full_df, "real_data.parquet")
full_df <- arrow::read_parquet("real_data.parquet")

exp_design <- data.frame(
    condition=c("FINE2a", "FINE2b", "TesR7A", "TesR7B"),
    DIUID = c("FINE", "FINE", "TesR", "TesR")
)
reference_transcripts_data <- readr::read_tsv(
    "Homo_sapiens.GRCh38.105.gtf.transcripts.tsv",
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
)

# cnames <- grep("dge1", full_df %>% colnames(), value = TRUE)
cnames <- full_df %>% colnames()
cnames1 <- grep("FINE", cnames, value = TRUE)
cnames2 <- grep("TesR", cnames, value = TRUE)

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
clusterExport(cl, varlist = ls())

retl <- parSapply(
        cl,
        unique(gene_transcript_df_detected$GENE_ID),
        function(gene_id){
            library(tidyverse)
            this_gene_transcript_table <- list()
            this_transcript_ids <- gene_transcript_df_detected %>%
                dplyr::filter(GENE_ID == gene_id) %>%
                dplyr::select(TRANSCRIPT_ID) %>%
                unlist() %>%
                as.character()
            if (length(this_transcript_ids) >= 2){
                this_gene_transcript_table[[gene_id]] <- this_transcript_ids
            }
            rm(gene_id, this_transcript_ids)
            this_gene_transcript_table
        }
    )
gene_transcript_table <- Reduce(c, retl)

rm(gene_transcript_df_detected)

full_df_mutated <- as.data.frame(full_df)
rownames(full_df_mutated) <- full_df_mutated$TRANSCRIPT_ID
full_df_mutated$TRANSCRIPT_ID <- NULL

c_df <- data.frame()
gene_ids <- names(gene_transcript_table)


for (i in seq_along(gene_ids)) {
    gene_id <- gene_ids[i]
    print(sprintf("%s -- %d/%d", gene_id, i, length(gene_ids)))
    sel_df <- full_df_mutated[gene_transcript_table[[gene_id]], ]
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
    )) +
    geom_label(aes(
        x=l2fc,
        y=basemean,
       label=TRANSCRIPT_ID
    )) +
    scale_color_manual(
        values = c("black", "red")
    ) +
    theme_bw()

ggplot(c_df) +
    geom_point(aes(
        x=l2fc,
        y=log10(pv),
        color=pv < 0.05 & abs(l2fc) > 2
    )) +
    scale_color_manual(
        values = c("black", "red")
    ) +
    theme_bw()
