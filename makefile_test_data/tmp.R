library(tidyverse)

all_table <- read_tsv("makefile_test_data/aaaaaaaaa.tsv")

fig_data <- all_table %>%
     dplyr::select(YASIM_1_ACTUAL_RPKM, FEATURECOUNTS_1_ACTUAL_RPKM, CPPTETGS_1_ACTUAL_RPKM) %>%
     tidyr::gather(key="key", value = "value", -YASIM_1_ACTUAL_RPKM)

g <- ggplot(fig_data)+
    geom_point(
        aes(x=YASIM_1_ACTUAL_RPKM, y=value, color=key),
        alpha=0.5
    )
ggsave("1.png",g, height = 10, width = 15)
