library("EBImage")
library("MASS")


sim_grid <- function(dims, grid_lines=10, plot=FALSE) {
    w <- dims[1]
    h <- dims[2]

    out <- matrix(rep(0, w*h), nrow=h)

    .w <- rep(floor(w/grid_lines), grid_lines) * c(1:grid_lines)
    .h <- rep(floor(h/grid_lines), grid_lines) * c(1:grid_lines)


    .x <- sapply(.w, function(x) {
        out[, x] <<- rep(1, nrow(out))
    })

    .x <- sapply(.h, function(x) {
        out[x, ] <<- rep(1, ncol(out))
    })

    out <- EBImage::Image(t(out), dim=dims, colormode="gray")

    if (plot) {
        quartz()
        display(out, method="raster", interpolate=FALSE)
    }

    out
}

out <- sim_grid(c(10, 10), 5, FALSE)
display(out, method="raster", interpolate=FALSE)


sim_transform <- function(img, A, t) {
    w <- dim(img)[1]
    h <- dim(img)[2]

    img <- t(img@.Data)

    Ainv <- MASS::ginv(A)
    trial <- matrix(c(1,1, 1,w, h,1, h,w), nrow=4, byrow=TRUE)

    y <- c()

    .x <- apply(trial, 1, function(x) {
        y <<- c(y, Ainv %*% (matrix(x, nrow=2) - t))
    })

    y <- matrix(y, ncol=2, byrow=TRUE)

    xmin <- min(1, floor(min(y[, 1])))
    xmax <- max(w, ceiling(max(y[, 1])))
    ymin <- min(1, floor(min(y[, 2])))
    ymax <- max(h, ceiling(max(y[, 2])))

    out <- sapply(xmin:xmax, function(x, ...) {
        sapply(ymin:ymax, function(y, ...) {
            x <- ...[[1]]
            img <- ...[[2]]
            A <- ...[[3]]
            t <- ...[[4]]

            .x <- A %*% matrix(c(x,y), nrow=2) + t

            if (.x[1] <=0 | .x[1] > h | .x[2] <= 0 | .x[2] > w) {
                # c(0, 0)
                0
            } else {
                # .x
                img[.x[1], .x[2]]
            }

        }, c(list(x), ...))
    }, list(img, A, t))

    out <- as.Image(t(out))
    out
}


out <- sim_grid(c(100, 100), 5, FALSE)
t <- matrix(c(0, 0), nrow=2)
A <- matrix(c(1.3, 0.2, 0.2, 0.7), nrow=2, byrow=TRUE)
tout <- sim_transform(out, A, t)

display(out, "raster")
display(tout, "raster")
