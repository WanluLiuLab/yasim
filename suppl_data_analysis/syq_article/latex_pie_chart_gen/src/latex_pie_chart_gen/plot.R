library(argparser)

p <- argparser::arg_parser("")
p <- argparser::add_argument(p, "--src_data_csv_file_path", "")
p <- argparser::add_argument(p, "--dst_fig_dir_path", "")
argv <- argparser::parse_args(p)

library(tidyverse)
library(ggpubr)
library(parallel)
cl <- parallel::makeCluster(spec = 8)

dir.create(argv$dst_fig_dir_path, showWarnings = FALSE)

data <- readr::read_csv(
    argv$src_data_csv_file_path,
    col_types = c(
        Software = col_character(),
        Dataset = col_character(),
        .default = col_double()
    )
)

data <- data %>%
    dplyr::mutate(
        fig_filename = file.path(
            argv$dst_fig_dir_path,
            sprintf("%s-%s.pdf", Software, Dataset)
        )
    )

data_long <- data %>%
    tidyr::gather(
        key = "type",
        value = "count",
        -Software,
        -Dataset,
        -fig_filename
    )

clusterExport(cl, varlist = ls(), envir = environment())
parLapply(cl, unique(data_long$fig_filename), function(fn) {
    colors <- c("#B2182B", "#F4A582", "#9E9AC8", "#92C5DE", "#2166AC")
    g <- ggpubr::ggpie(
        dplyr::filter(
            data_long,
            fig_filename == fn
        ),
        "count",
        label = "type",
        fill = "type",
        color = "white",
    ) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_fill_manual(values = colors)
    ggplot2::ggsave(fn, g, width = 5, height = 5)
})

stopCluster(cl)
