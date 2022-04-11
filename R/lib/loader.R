load_package <- function(name) {
    message(sprintf("Loading %s", name))
    suppressWarnings(suppressMessages(library(name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)))
}

load_package("tidyverse")
