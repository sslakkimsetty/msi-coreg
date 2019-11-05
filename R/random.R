

construct_sample <- function(dims, size=50) {
    w <- dims[1]
    h <- dims[2]

    out <- matrix(
        c(sample.int(w, size=size), sample.int(n=h, size=size)),
        nrow = size
    )

    out
}

