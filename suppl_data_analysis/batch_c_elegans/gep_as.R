library("tidyverse")
library("svglite")

as_event_types <- c("alt_3prime", "alt_5prime", "exon_skip", "intron_retention", "mult_exon_skip", "mutex_exons")

#' Get dataset metadata
datasets_metadata <-
    readr::read_csv("datasets.csv", show_col_types = FALSE) %>%
        dplyr::filter(! Read %in% c("SRR8081024", "SRR8081025"), Generation != "TGS") %>%
        dplyr::mutate(filename = sprintf("./datasets/%s/%s_SPLADDER/all_as.tsv.aggregated.tsv", Model, Read)) %>%
        dplyr::mutate(
            alt_3prime=c(), alt_5prime=c(), exon_skip=c(), intron_retention=c(), mult_exon_skip=c(), mutex_exons=c()
        )

#' Raw sequence run accessions
all_run_accessions <- unlist(datasets_metadata$Read)

#' Read data
for (i in seq_len(nrow(datasets_metadata))) {
    data_tuple <- datasets_metadata[i,]
    if (file.exists(data_tuple$filename)){
        x <- readr::read_tsv(data_tuple$filename, show_col_types = FALSE)
        for (as_event_type in as_event_types){
            datasets_metadata[[as_event_type]][i]<-sum(x[[as_event_type]])
        }
        rm(x)
    } else{
        warning(sprintf("%s not found!", data_tuple$filename))
    }
    rm(data_tuple)
}

datasets_metadata_long <- datasets_metadata %>%
    dplyr::select(!(Generation:Project)) %>%
    dplyr::select(!(filename)) %>%
    tidyr::gather(key = "ASEventType", value = "Number", -Read)

g <- ggplot(datasets_metadata_long) +
    geom_bar(aes(x=Read, y=Number, fill=ASEventType), stat="identity", position = "fill") +
    ggtitle("Ratio between alternative splicing events") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave(
    "ratio_as_events.svg",
    g,
    width=8,height=5
)

datasets_metadata_sum <- datasets_metadata %>%
    dplyr::select(all_of(as_event_types)) %>%
    summarise_all(sum)
