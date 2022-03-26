library(tidyverse)
library(pheatmap)
library(ComplexUpset)

all_as <- read_tsv("all_as.tsv", show_col_types = FALSE, col_names = c("GENE_ID", "EVENT"))

all_as_summary <- all_as %>%
    dplyr::count(GENE_ID, EVENT)

events <- unique(all_as_summary$EVENT)

all_as_summary_wide <- tidyr::spread(
    all_as_summary, key=EVENT, value=n, fill=0
) %>%
    dplyr::rowwise()%>%
    dplyr::mutate(sum=sum(c_across(cols = events)))

upset(data = all_as_summary_wide, intersect = events)

all_as_summary_wide_pheatmap <- all_as_summary_wide %>%
    dplyr::arrange(desc(sum)) %>%
    head(n = 500) %>%
    dplyr::select(!sum)

gene_id_t500 <- all_as_summary_wide_pheatmap$GENE_ID

all_as_summary_wide_pheatmap_matrix <- all_as_summary_wide_pheatmap %>%
    dplyr::select(!GENE_ID) %>%
    as.matrix()

rownames(all_as_summary_wide_pheatmap_matrix) <- gene_id_t500
