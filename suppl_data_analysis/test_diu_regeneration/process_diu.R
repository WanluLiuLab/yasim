library(argparser)
library(tidyverse)
library(IsoformSwitchAnalyzeR)

library(parallel)


read_fc <- function(fname){
    readr::read_tsv(
        fname,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            Chr = col_character(),
            Start = col_character(),
            End = col_character(),
            Strand = col_character(),
            Length = col_integer(),
            NumReads = col_integer()
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
}

fine2a <- read_fc(
    "FINE2a.fastq.gz.sam.bam.chr21.bam.fC.txt"
) %>%
    dplyr::rename(
        FINE2a=NumReads
    )

fine2b <- read_fc(
    "FINE2b.fastq.gz.sam.bam.chr21.bam.fC.txt"
) %>%
    dplyr::rename(
        FINE2b=NumReads
    )

tesr7a <- read_fc(
    "TesR7A.fastq.gz.sam.bam.chr21.bam.fC.txt"
) %>%
    dplyr::rename(
        TesR7A=NumReads
    )

tesr7b <- read_fc(
    "TesR7B.fastq.gz.sam.bam.chr21.bam.fC.txt"
) %>%
    dplyr::rename(
        TesR7B=NumReads
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

exp_design <- data.frame(
    sampleID=c("FINE2a", "FINE2b", "TesR7A", "TesR7B"),
    condition = c("FINE", "FINE", "TesR", "TesR")
)

full_df_isar <- full_df %>% as.data.frame() %>% dplyr::rename(isoform_id=TRANSCRIPT_ID)


### Create switchAnalyzeRlist
aSwitchList <- IsoformSwitchAnalyzeR::importRdata(
    isoformNtFasta = "chr21_trans.fa",
    showProgress = TRUE,
    isoformCountMatrix = full_df_isar,
    designMatrix = exp_design,
    isoformExonAnnoation = "hg38.ncbiRefSeq.chr21.gtf",
    ignoreAfterPeriod = FALSE
)

SwitchListAnalyzed <- IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchList,
    reduceToSwitchingGenes = TRUE,
    reduceFurtherToGenesWithConsequencePotential = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE
) %>%
    IsoformSwitchAnalyzeR::analyzeORF(showProgress = FALSE) %>%
    IsoformSwitchAnalyzeR::extractSequence() %>%
    IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(showProgress = FALSE, onlySwitchingGenes = FALSE)


exampleSwitchListAnalyzed <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
    SwitchListAnalyzed,
    consequencesToAnalyze = c('intron_retention', 'NMD_status', 'ORF_seq_similarity'),
    dIFcutoff = 0.1,
    alpha = 0.05,
    showProgress = TRUE
)

##Obtain all switching genes
All_significant_DIU <- IsoformSwitchAnalyzeR::extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = NA,
    extractGenes = FALSE,
    sortByQvals = TRUE
)

