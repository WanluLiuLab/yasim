file_description <- ""

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--ss_tsv", help = ".tsv file produced by get_sam_satistics.py", type = "character")
p <- add_argument(p, "--output", help = "Output basename.", type = "character")
argv <- parse_args(p)

source(argv$libfile)
ss_tsv <- argv$ss_tsv

ss_data <- read_tsv(ss_tsv, col_types = ss_tsv_col_types)


g <- ggplot(ss_data) + geom_density(aes(x = LEN))
ggsave(paste(argv$output, "LEN", "png", sep = "."), plot = g)

g <- ggplot(ss_data) + geom_density(aes(x = GC))
ggsave(paste(argv$output, "GC", "png", sep = "."), plot = g)

# if (1 == length(unique(ss_data$LEN))){
    g <- ggplot(ss_data) + geom_density(aes(x = MAPQ))
    ggsave(paste(argv$output, "MAPQ", "png", sep = "."), plot = g)

    g <- ggplot(ss_data) + geom_density(aes(x = ALNQ))
    ggsave(paste(argv$output, "ALNQ", "png", sep = "."), plot = g)
# } else {
#     g <- ggplot(ss_data) + geom_density2d_filled(aes(x = LEN, y=MAPQ))
#     ggsave(paste(argv$output, "LEN-MAPQ", "png", sep = "."), plot = g)
#
#     g <- ggplot(ss_data) + geom_density2d_filled(aes(x = LEN, y=ALNQ))
#     ggsave(paste(argv$output, "LEN-ALNQ", "png", sep = "."), plot = g)
# }

# g <- ggplot(ss_data) + geom_density2d_filled(aes(x = GC, y=MAPQ))
# ggsave(paste(argv$output, "GC-MAPQ", "png", sep = "."), plot = g)
#
# g <- ggplot(ss_data) + geom_density2d_filled(aes(x = GC, y=ALNQ))
# ggsave(paste(argv$output, "GC-ALNQ", "png", sep = "."), plot = g)
