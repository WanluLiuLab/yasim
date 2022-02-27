file_description <- ""

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--htseq_quant_tsv", help = ".tsv file produced by htseq-count.", type = "character")
p <- add_argument(p, "--yasim_tsv", help = "Gene depth .tsv file produced by tasim.", type = "character")
p <- add_argument(p, "--output", help = "Output basename.", type = "character")
argv <- parse_args(p)

source(argv$libfile)
htseq_quant_tsv <- argv$htseq_quant_tsv
yasim_tsv <- argv$yasim_tsv

yasim_data <- read_tsv(yasim_tsv, col_types = yasim_depth_tsv_col_types)
htseq_quant_data <- read_tsv(
    htseq_quant_tsv,
    col_types = htseq_quant_tsv_col_types,
    col_names = htseq_quant_tsv_col_names
) %>% dplyr::filter(NumReads > 0)

all_table <- dplyr::inner_join(yasim_data, htseq_quant_data, by = c("TRANSCRIPT_ID" = "Name"))

message(sprintf("Read %d from yasim_tsv and %d from htseq_quant_tsv. %d left merged.",
                nrow(yasim_data), nrow(htseq_quant_data), nrow(all_table)))

g <- ggplot(all_table) + stat_summary(aes(x = DEPTH, y = NumReads))
ggsave(paste(argv$output, "png", sep = "."), plot = g)
