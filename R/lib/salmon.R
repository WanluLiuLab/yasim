salmon_quant_sf_col_types <- cols(
    Name = col_character(),
    Length = col_double(),
    EffectiveLength = col_double(),
    TPM = col_double(),
    NumReads = col_double()
)

get_salmon_data <- function(salmon_quant_sf, n, ATTR_NAME) {
    message(sprintf("Reading %s...", salmon_quant_sf))
    salmon_quant_sf_data <- read_tsv(salmon_quant_sf, col_types = salmon_quant_sf_col_types) %>%
        dplyr::transmute(
            !!rlang::sym(ATTR_NAME) := Name,
            !!rlang::sym(sprintf("SALMON_%d_ACTUAL_N_OF_READS", n)) := NumReads
        )
    message(sprintf("Reading %s... DONE", salmon_quant_sf))
    return(salmon_quant_sf_data)
}
