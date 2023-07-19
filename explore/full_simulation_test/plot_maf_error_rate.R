library("tidyverse")

all_error <- readr::read_tsv(
    "all_last_mapq.tsv",
    col_types = c(
        FILENAME = col_character(),
        INSERTION = col_double(),
        DELETION = col_double(),
        MATCH = col_double(),
        SUBSTITUTION = col_double()
    )
) %>%
    tidyr::gather(
        key = "EventType",
        value = "EventCount",
        -FILENAME
    ) %>%
    dplyr::mutate(
        Condition = FILENAME %>%
            stringr::str_replace(
                "ce11_",
                ""
            ) %>%
            stringr::str_replace(
                ".maf.gz",
                ""
            )
    ) %>%
    tidyr::separate(
        "Condition",
        c("SIMULATOR", "SEQUENCER", "MODE", "DGEID", "DIUID", "REPID")
    )

g <- ggplot(all_error) +
    geom_bar(
        aes(
            y = sprintf("%s_%s_%s", SIMULATOR, SEQUENCER, MODE),
            x = EventCount,
            fill = EventType
        ),
        position = "fill",
        stat = "identity"
    ) +
    facet_wrap(DGEID ~ DIUID) +
    theme_bw()
ggsave("maf_error_rate.pdf", g, width = 8, height = 5)

