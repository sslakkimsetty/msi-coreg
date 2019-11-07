

estimate_gaussian <- function(samp_int) {
    mu <- mean(samp_int)
    var <- var(samp_int)

    list(mu=mu, var=var)
}


g_density <- function(zi, zj, gaussian_params) {
    sigma <- gaussian_params[[2]]

    Gphi <- (2*pi)^(-1/2) * sigma^(-1/2)
    Gphi <- Gphi * exp((-1/2) * (zi-zj) * (1/sigma) * (zi-zj))

    Gphi
}


g_weighting_factor <- function(sampleA, zi, zj, gaussian_params) {
    # If sampleA is a (x-y-intensities) matrix, extract intensities
    if (is.matrix(sampleA)) {
        sampleA <- sampleA[, 3]
    }

    num <- g_density(zi, zj, gaussian_params)
    den <- sapply(sampleA, function(x, ...) {
        zi <- ...[[1]]
        gaussian_params <- ...[[2]]
        g_density(zi, x, gaussian_params)},
        list(zi, gaussian_params))

    out <- (num / sum(den))

    out
}


