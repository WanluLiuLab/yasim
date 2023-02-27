library("tidyverse")

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
        SEQUENCER = condition_break[[1]][2],
        MODE = condition_break[[1]][3],
        DGEID = condition_break[[1]][4],
        DIUID = condition_break[[1]][5],
        REPID = condition_break[[1]][6],
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
save(
    df_list = df_list,
    file = "dge_df_list_prepared.rd.xz",
    compress = "xz",
    compression_level = 9
)
