file_description <- "summarize_as.R -- Suumarize alternative splicing information by Spladder"

library(argparser)

p <- arg_parser(file_description)
p <- add_argument(p, "--libfile", help = "Libfile to read.", type = "character")
p <- add_argument(p, "--filename", help = "all_as.tsv produced by extract_as_from_spladder.py", type = "character")
argv <- parse_args(p)

source(argv$libfile)

load_package("pheatmap")
load_package("ComplexUpset")

filename <- argv$filename

all_as <- read_tsv(filename, show_col_types = FALSE, col_names = c("GENE_ID", "EVENT"))
all_as_summary <- all_as %>%
    dplyr::count(GENE_ID, EVENT)

events <- unique(all_as_summary$EVENT)

all_as_summary_wide <- tidyr::spread(
    all_as_summary, key = EVENT, value = n, fill = 0
) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sum = sum(c_across(cols = all_of(events))))

write_tsv(all_as_summary_wide, paste0(filename, ".aggregated.tsv"))

all_as_summary_wide_upset <- all_as_summary_wide
all_as_summary_wide_upset[events] <- all_as_summary_wide[events] != 0

p <- upset(data = all_as_summary_wide_upset, intersect = events)
ggsave(
    paste0(filename, ".upset.pdf"),
    p,
    width = 16,
    height = 10
)

all_as_summary_wide_pheatmap <- all_as_summary_wide %>%
    dplyr::arrange(desc(sum)) %>%
    dplyr::select(!sum)

gene_id <- all_as_summary_wide_pheatmap$GENE_ID

all_as_summary_wide_pheatmap_matrix <- all_as_summary_wide_pheatmap %>%
    dplyr::select(!GENE_ID) %>%
    as.matrix()

rownames(all_as_summary_wide_pheatmap_matrix) <- gene_id

pdf(
    paste0(filename, ".top50_heatmap.pdf"),
    width = 10, height = 10
)
pheatmap(
    head(all_as_summary_wide_pheatmap_matrix, 50),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    angle_col = 0
)
while (!is.null(dev.list())) dev.off()
