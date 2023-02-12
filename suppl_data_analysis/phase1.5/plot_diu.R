# Import the required packages
library(IsoformSwitchAnalyzeR)
library(tidyverse)

all_fc_data_isoform_level <- NULL
experiment_design <- NULL
fns <- Sys.glob("ce11_*.fq.gz.bam.fc_stringtie.tsv")
conditions <- fns %>%
    stringr::str_replace(".fq.gz.bam.fc_stringtie.tsv", "") %>%
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
    if (is.null(all_fc_data_isoform_level)) {
        all_fc_data_isoform_level <- this_fc_data_isoform_level %>%
            dplyr::select(TRANSCRIPT_ID)
    }
    all_fc_data_isoform_level <- all_fc_data_isoform_level %>%
        dplyr::inner_join(this_fc_data_isoform_level, by = "TRANSCRIPT_ID")

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
        fc_data_fn,
        this_experiment_design,
        i
    )
    gc()
}

arrow::write_parquet(all_fc_data_isoform_level, "diu_fc_data.parquet")

isar_design_matrix <- experiment_design %>%
    dplyr::filter(SIMULATOR=="pbsim3") %>%
    dplyr::filter(DIUID=="diu1") %>%
    dplyr::filter(REPID=="0") %>%
    dplyr::transmute(sampleID=condition, condition=DGEID)
isar_count_matrix <- all_fc_data_isoform_level %>%
    dplyr::select(tidyselect::contains("pbsim3")) %>%
    dplyr::select(tidyselect::contains("diu1")) %>%
    dplyr::select(tidyselect::ends_with("0")) %>%
    dplyr::mutate(isoform_id=all_fc_data_isoform_level$TRANSCRIPT_ID) %>%
    dplyr::relocate(isoform_id, .before = tidyselect::where(is.numeric))

### Create switchAnalyzeRlist
aSwitchList <- IsoformSwitchAnalyzeR::importRdata(
    isoformNtFasta = "ce11_trans_stringtie.fa",
    showProgress = TRUE,
    isoformCountMatrix = isar_count_matrix,
    designMatrix = isar_design_matrix,
    isoformExonAnnoation = "stringtie_merged.gtf",
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
    IsoformSwitchAnalyzeR::analyzeORF(showProgress = TRUE) %>%
    IsoformSwitchAnalyzeR::extractSequence() %>%
    IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(showProgress = TRUE, onlySwitchingGenes = FALSE)


exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
    SwitchListAnalyzed,
    consequencesToAnalyze = c('intron_retention', 'NMD_status', 'ORF_seq_similarity'),
    dIFcutoff = 0.1,
    alpha = 0.05,
    showProgress = TRUE
)

##Obtain all switching genes
All_significant_DIU <- extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = NA,
    extractGenes = FALSE,
    sortByQvals = TRUE
)
