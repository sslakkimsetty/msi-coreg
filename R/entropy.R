

estimate_gaussian <- function(sample) {
    mu <- mean(sample[, 3])
    var <- var(sample[, 3])

    list(mu=mu, var=var)
}


superpose_gaussian <- function(sampleA)


