library("tidyverse")
library("pheatmap")
library("scales")
library("ggridges")
library("arrow")
library("corrplot")
library("DESeq2")

all_fc_data <- NULL
experiment_design <- NULL
for (fc_data_fn in Sys.glob("ce11_*.fq.bam.fc.tsv")) {
    condition <- fc_data_fn %>%
        stringr::str_replace(".fq.bam.fc.tsv", "") %>%
        stringr::str_replace("ce11_", "")
    this_fc_data <- readr::read_tsv(
        fc_data_fn,
        col_types = cols(
            TRANSCRIPT_ID = col_character(),
            Chr = col_character(),
            Start = col_character(),
            End = col_character(),
            Strand = col_character(),
            Length = col_number(),
            NumReads = col_double(),
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
    if (is.null(all_fc_data)) {
        all_fc_data <- this_fc_data %>%
            dplyr::select(TRANSCRIPT_ID)
    }
    this_fc_data <- this_fc_data %>%
        dplyr::transmute(
            TRANSCRIPT_ID = TRANSCRIPT_ID,
            !!rlang::sym(condition) := NumReads
        )
    all_fc_data <- all_fc_data %>%
        dplyr::inner_join(this_fc_data, by = "TRANSCRIPT_ID")

    condition_break <- strsplit(condition, "_")
    this_experiment_design <- data.frame(
        model=condition_break[[1]][1],
        group=condition_break[[1]][3],
        condition=condition
    )
    if (is.null(experiment_design)){
        experiment_design <- this_experiment_design
    } else{
        experiment_design <- experiment_design %>%
            dplyr::rows_append(
            this_experiment_design
        )
    }
    rm(condition, this_fc_data, fc_data_fn, this_experiment_design)
}

experiment_design <- as.data.frame(experiment_design)
rownames(experiment_design) <- experiment_design$condition
experiment_design$condition <- NULL

all_fc_data <- as.data.frame(all_fc_data)
rownames(all_fc_data) <- all_fc_data$TRANSCRIPT_ID
all_fc_data$TRANSCRIPT_ID <- NULL

dds <- DESeqDataSetFromMatrix(
    countData = all_fc_data,
    colData = experiment_design,
    design = ~ model + group
)
dds <- DESeq(dds)
