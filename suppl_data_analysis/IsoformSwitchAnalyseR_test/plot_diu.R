# Import the required packages
library(IsoformSwitchAnalyzeR)
library(tidyverse)

isar_count_matrix <- data.frame(
    S1.1=c(1, 4, 4, 70),
    S2.1=c(1, 4, 70, 4),
    S1.2=c(2, 3, 5, 70),
    S2.2=c(2, 3, 71, 5)
)
row.names(isar_count_matrix) <- c("F1.2", "F1.1", "F2.2", "F2.1")
isar_design_matrix <- data.frame(
    sampleID=c("S1.1", "S2.1", "S1.2", "S2.2"),
    condition=c("S1", "S2", "S1", "S2")
)


### Create switchAnalyzeRlist
aSwitchList <- IsoformSwitchAnalyzeR::importRdata(
    isoformNtFasta = "fake_trans.fa",
    showProgress = TRUE,
    isoformCountMatrix = isar_count_matrix,
    designMatrix = isar_design_matrix,
    isoformExonAnnoation = "fake.gtf",
    ignoreAfterPeriod = FALSE #if using the reference gtf, then this should set TRUE.
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
