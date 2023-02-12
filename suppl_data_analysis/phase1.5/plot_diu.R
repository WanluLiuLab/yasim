# Import the required packages
library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(arrow)

all_fc_data_isoform_level <- arrow::read_parquet("diu_fc_data_gt.parquet")
experiment_design <- arrow::read_parquet("dge_experiment_design.parquet")

isar_design_matrix <- experiment_design %>%
    dplyr::filter(SIMULATOR=="pbsim3") %>%
    dplyr::filter(REPID=="0") %>%
    dplyr::transmute(sampleID=condition, condition=DGEID)
isar_count_matrix <- all_fc_data_isoform_level %>%
    dplyr::select(tidyselect::contains("pbsim3")) %>%
    dplyr::select(tidyselect::ends_with("0")) %>%
    dplyr::mutate(isoform_id=all_fc_data_isoform_level$TRANSCRIPT_ID) %>%
    dplyr::relocate(isoform_id, .before = is.numeric)

### Create switchAnalyzeRlist
aSwitchList <- IsoformSwitchAnalyzeR::importRdata(
    isoformNtFasta = "ce11_trans_as.chr1.fa",
    showProgress = TRUE,
    isoformCountMatrix = isar_count_matrix,
    designMatrix = isar_design_matrix,
    isoformExonAnnoation = "ce11.ncbiRefSeq_as.chr1.gtf",
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
