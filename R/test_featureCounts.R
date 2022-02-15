file_description <- ""

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--featureCounts_tsv", help = ".tsv file produced by featureCounts", type = "character")
p <- add_argument(p, "--yasim_tsv", help = "Gene depth .tsv file produced by tasim.", type = "character")
p <- add_argument(p, "--output", help = "Output basename", type = "character")
argv <- parse_args(p)

source(argv$libfile)
featureCounts_tsv <- argv$featureCounts_tsv
yasim_tsv <- argv$yasim_tsv

yasim_data <- read_tsv(yasim_tsv, , col_types = yasim_tsv_col_types)
featureCounts_data <- read_tsv(
    featureCounts_tsv,
    quote = "\'",
    col_types = featureCounts_tsv_col_types,
    col_names = featureCounts_tsv_col_names,
    comment="#"
) %>% dplyr::filter(NumReads>0)
print(head(featureCounts_data))
all_table <-  dplyr::inner_join(yasim_data, featureCounts_data, by=c("gene_name" = "Geneid"))

message(sprintf("Read %d from yasim_tsv and %d from salmon_quant_sf. %d left merged.",
                nrow(yasim_data), nrow(featureCounts_data), nrow(all_table)))

g <- ggplot(all_table) + stat_summary(aes(x=depth,y=NumReads))
ggsave(paste0(argv$output,".png"), plot=g)
