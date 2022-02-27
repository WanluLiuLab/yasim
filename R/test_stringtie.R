file_description <- ""

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--stringtie_quant_tsv", help = ".tsv file produced by parse_stringtie_into_tsv.py", type = "character")
p <- add_argument(p, "--yasim_tsv", help = "Gene depth .tsv file produced by tasim.", type = "character")
p <- add_argument(p, "-e", help = "Whether stringtie was called with `-e` option.", flag = TRUE)
p <- add_argument(p, "--output", help = "Output basename", type = "character")
argv <- parse_args(p)

source(argv$libfile)
stringtie_quant_tsv <- argv$stringtie_quant_tsv
yasim_tsv <- argv$yasim_tsv

yasim_data <- read_tsv(yasim_tsv, , col_types = yasim_depth_tsv_col_types)
stringtie_quant_data <- read_tsv(stringtie_quant_tsv, quote = "\'", col_types = stringtie_quant_tsv_col_types)

if (argv$e) {
    get_id <- "transcript_id"
} else {
    get_id <- "reference_id"
}

all_table <- dplyr::inner_join(yasim_data, stringtie_quant_data, by = c("TRANSCRIPT_ID" = get_id))

message(sprintf("Read %d from yasim_tsv and %d from salmon_quant_sf. %d left merged.",
                nrow(yasim_data), nrow(stringtie_quant_data), nrow(all_table)))

g <- ggplot(all_table) + stat_summary(aes(x = DEPTH, y = cov))
ggsave(paste0(argv$output, ".png"), plot = g)
