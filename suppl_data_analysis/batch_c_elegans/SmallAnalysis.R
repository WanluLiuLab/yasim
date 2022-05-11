library("tidyverse")
library("svglite")

#' Get dataset metadata
datasets_metadata <-
    readr::read_csv("datasets.csv", show_col_types = FALSE) %>%
        dplyr::mutate(filename = sprintf("./datasets/%s/%s.transcript.depth.tsv", Model, Read))

#' Raw sequence run accessions
all_run_accessions <- unlist(datasets_metadata$Read)

#' All GEP data in depth
all_gep_data <- NA

#' Read data
for (i in seq_len(nrow(datasets_metadata))) {
    data_tuple <- datasets_metadata[i,]
    if (!tibble::is_tibble(all_gep_data)) {
        all_gep_data <-
            readr::read_tsv(data_tuple$filename, show_col_types = FALSE) %>%
                dplyr::mutate(!!rlang::sym(data_tuple$Read) := AVG_DEPTH) %>%
                dplyr::select(!AVG_DEPTH)
    } else {
        all_gep_data <-
            readr::read_tsv(data_tuple$filename, show_col_types = FALSE) %>%
                dplyr::select(TRANSCRIPT_ID, AVG_DEPTH) %>%
                dplyr::filter(AVG_DEPTH > 0) %>%
                dplyr::inner_join(all_gep_data, by = "TRANSCRIPT_ID") %>%
                dplyr::mutate(!!rlang::sym(data_tuple$Read) := AVG_DEPTH) %>%
                dplyr::select(!AVG_DEPTH)
    }
    rm(data_tuple)
}

#' Get number of records
n_of_records <- nrow(all_gep_data)

#' Normalize data by scaling them into (-1, 1)
all_gep_data[, 9:length(all_gep_data)] <-
    scale(all_gep_data[, 9:length(all_gep_data)])

#' Transform into long form
all_gep_data_long <- all_gep_data %>%
    dplyr::select(!(GENE_ID:GC)) %>%
    tidyr::gather(key = "RunAccession", value = "Depth", -TRANSCRIPT_ID)


g <- ggplot(all_gep_data_long) +
    geom_density(aes(x = Depth, color = RunAccession)) +
    scale_x_continuous(name = "log(NormalizedDepth)", trans = "log") +
    scale_y_continuous(name = "density", limits = c(0, 0.35)) +
    theme_bw() +
    ggtitle("Normalized depth density distribution accross all samples")
ggsave("gep_density.svg", g, width = 8, height = 5)

#' Fit Zipf's law
#' Make data headless, which is also used for aggregation
#' and sorting.
all_gep_data_zipf_depth <- all_gep_data %>%
    dplyr::select(!(TRANSCRIPT_ID:GC))

#' Sort lines
for (i in all_run_accessions) {
    all_gep_data_zipf_depth[[i]] <- sort(all_gep_data_zipf_depth[[i]])
}

#' Transform to long form
#' while decreasing data size
all_gep_data_zipf_depth_long <- all_gep_data_zipf_depth %>%
    dplyr::mutate(RANK = n_of_records - 1:n_of_records) %>%
    tidyr::gather(key = "RunAccession", value = "Depth", -RANK)

g <- ggplot(all_gep_data_zipf_depth_long) +
    geom_point(aes(y = Depth, x = RANK, color = RunAccession), size=0.2) +
    scale_x_continuous(name = "log(Rank)", trans = "log") +
    scale_y_continuous(name = "log(NormalizedDepth)", trans = "log") +
    theme_bw() +
    ggtitle("Normalized depth for Zipf's Law accross all samples")
ggsave("gep_zipf.svg", g, width = 8, height = 5)

#' Aggregate by sequencers
all_gep_data_sequencer <- all_gep_data_zipf_depth %>%
    dplyr::mutate(RANK = 1:dim(all_gep_data_zipf_depth)[1])

for (model_name in unique(datasets_metadata$Model)) {
    reads <- datasets_metadata %>%
        dplyr::filter(Model == model_name) %>%
        dplyr::select(Read) %>%
        unlist()
    all_gep_data_sequencer <- all_gep_data_sequencer %>%
        dplyr::mutate(!!model_name := rowMeans(select(., all_of(reads)))) %>%
        dplyr::select(!all_of(reads))
}

all_gep_data_sequencer_long <- all_gep_data_sequencer %>%
    tidyr::gather(key = "Sequencer", value = "Depth", -RANK)

g <- ggplot(all_gep_data_sequencer_long) +
    geom_density(aes(x = Depth, color = Sequencer)) +
    scale_x_continuous(name = "log(NormalizedDepth)", trans = "log") +
    scale_y_continuous(name = "density", limits = c(0, 0.35)) +
    theme_bw() +
    ggtitle("Normalized depth density distribution accross all sequencers")
ggsave("gep_density_sequencer.svg", g, width = 8, height = 5)


#' Aggregate by generation
all_gep_data_generation <- all_gep_data_zipf_depth %>%
    dplyr::mutate(RANK = 1:dim(all_gep_data_zipf_depth)[1])


for (generation_name in unique(datasets_metadata$Generation)) {
    reads <- datasets_metadata %>%
        dplyr::filter(Generation == generation_name) %>%
        dplyr::filter(Model != "NextSeq550") %>%
        dplyr::select(Read) %>%
        unlist()
    all_gep_data_generation <- all_gep_data_generation %>%
        dplyr::mutate(!!generation_name := rowMeans(select(., all_of(reads)))) %>%
        dplyr::select(!all_of(reads))
}

all_gep_data_generation_long <- all_gep_data_generation %>%
    dplyr::select(!dplyr::starts_with("SRR")) %>%
    tidyr::gather(key = "Generation", value = "Depth", -RANK)

g <- ggplot(all_gep_data_generation_long) +
    geom_density(aes(x = Depth, color = Generation)) +
    scale_x_continuous(name = "log(NormalizedDepth)", trans = "log") +
    scale_y_continuous(name = "density", limits = c(0, 0.35)) +
    theme_bw() +
    ggtitle("Normalized depth density distribution accross NGS and TGS")
ggsave("gep_density_generation.svg", g, width = 8, height = 5)
#
# q1 <- quantile(all_gep_data_generation$NGS, probs=seq(0, 1, 1/10E3))
# q2 <- quantile(all_gep_data_generation$TGS, probs=seq(0, 1, 1/10E3))
#
# g <- ggplot() +
#     geom_point(aes(x=q1, y=q2)) +
#     geom_abline(slope=1, intercept = 0) +
#     theme_bw() +
#     ggtitle("Quantile-Quantile Plot accross NGS and TGS")
# ggsave("gep_qq_generation.svg", g, width=8, height=5)
