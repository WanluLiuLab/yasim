library(tidyverse)
library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--stringtie_quant_tsv", help = ".tsv file produced by parse_stringtie_into_tsv.py", type = "character")
p <- add_argument(p, "--yasim_tsv", help = "Gene depth .tsv file produced by tasim.", type = "character")
p <- add_argument(p, "-e", help = "Whether stringtie was called with `-e` option.", type = "bool")
argv <- parse_args(p)

stringtie_quant_tsv <- argv$stringtie_quant_tsv
yasim_tsv <- argv$stringtie_quant_tsv

yasm <- read_tsv(yasim_tsv)
stringtie_quant <- read_tsv(stringtie_quant_tsv, quote = "\'")

if (argv$e){
    get_id <- "transcript_id"
} else{
    get_id <- "reference_id"
}


all_table <-  dplyr::inner_join(yasm, stringtie_quant, by=c("gene_name" = get_id))

g <- ggplot(all_table) + stat_summary(aes(x=depth,y=cov))
ggsave("yasim_to_stringtie_quant.png", plot=g)

