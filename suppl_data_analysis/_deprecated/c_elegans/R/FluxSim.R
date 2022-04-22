dFluxSim <- function(x, k, Y_max, a, b) {
    return(Y_max * x^k * exp(x / a * (x / b)^2))
}

pFluxSim <- function(q, k, Y_max, a, b) {
    if (length(q) == 0) {
        return(c())
    }
    return(sum(dFluxSim(c(1:q), k, Y_max, a, b)))
}


dZipf <- function(x, k, Y_max) {
    return(Y_max * x^k)
}

pZipf <- function(q, k, Y_max, a, b) {
    if (length(q) == 0) {
        return(NULL)
    }
    return(sum(dZipf(seq_along(q), k, Y_max)))
}

qqZipf <- function(x, k) {
    actual_data <- x
    Y_max <- max(x)
    theotical_data <- dZipf(seq_along(x), k, Y_max)
    qqplot(actual_data, theotical_data)
    abline(0, 1)
}

qqFluxSim <- function(x, k, a, b) {
    actual_data <- x
    Y_max <- max(x)
    theotical_data <- dFluxSim(seq_along(x), k, Y_max, a, b)
    qqplot(actual_data, theotical_data)
    abline(0, 1)
}
