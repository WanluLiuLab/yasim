library(tidyverse)
library(arrow)

all_fc_data_isoform_level <- NULL
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

    message(sprintf("Processing %d/%d", i, length(fns)))
    rm(
        condition,
        this_fc_data,
        this_fc_data_isoform_level,
        fc_data_fn,
        this_experiment_design,
        i
    )
    gc()
}

arrow::write_parquet(all_fc_data_isoform_level, "diu_fc_data.parquet")

all_fq_stats_isoform_level <- readr::read_tsv(
    "ce11.ncbiRefSeq_as.chr1.gtf.transcripts.tsv",
    col_types = c(
        TRANSCRIPT_ID = col_character(),
        GENE_ID = col_character(),
        NAIVE_LENGTH = col_integer(),
        TRANSCRIBED_LENGTH = col_integer(),
        EXON_NUMBER = col_integer()
    )
) %>%
    dplyr::select(TRANSCRIPT_ID)
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
    all_fq_stats_isoform_level <- all_fq_stats_isoform_level %>%
        dplyr::left_join(
            this_fq_stats_data_isoform_level,
            by = "TRANSCRIPT_ID"
        )

    message(sprintf("Processing %d/%d -- %s", i, length(fns), fns[i]))
    rm(
        this_fq_stats_data,
        this_fq_stats_data_isoform_level,
        condition,
        i
    )
    gc()
}

all_fq_stats_isoform_level <- all_fq_stats_isoform_level %>%
    replace(is.na(.), 0)

arrow::write_parquet(all_fq_stats_isoform_level, "diu_fq_stats.parquet")

all_fc_data_isoform_level <- NULL
fns <- Sys.glob("ce11_*.fq.gz.bam.fc_gt.tsv")
conditions <- fns %>%
    stringr::str_replace(".fq.gz.bam.fc_gt.tsv", "") %>%
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

    message(sprintf("Processing %d/%d", i, length(fns)))
    rm(
        condition,
        this_fc_data,
        this_fc_data_isoform_level,
        fc_data_fn,
        this_experiment_design,
        i
    )
    gc()
}

arrow::write_parquet(all_fc_data_isoform_level, "diu_fc_data_gt.parquet")

