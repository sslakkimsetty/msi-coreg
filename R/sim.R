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

    out <- EBImage::Image(t(out), dim = dims, colormode = "gray")

    if (plot) {
        quartz()
        display(out, method="raster", interpolate=FALSE)
    }

    out
}

out <- sim_grid(c(10, 10), 5, FALSE)
display(out, method="raster", interpolate=FALSE)


sim_affine <- function(img, Q) {
    w <- dim(img)[1]
    h <- dim(img)[2]

    img <- t(img@.Data)

    Ainv <- MASS::ginv(A)

    t <- matrix(c(t[2], t[1]), nrow=2)
    trial <- matrix(c(1,1, 1,w, h,1, h,w), nrow=4, byrow=TRUE)

    y <- c()

    .x <- apply(trial, 1, function(x) {
        y <<- c(y, Ainv %*% (matrix(x, nrow=2) - t))
    })

    y <- matrix(y, ncol=2, byrow=TRUE)

    xmin <- min(1, floor(min(y[, 1])))
    xmax <- max(h, ceiling(max(y[, 1])))
    ymin <- min(1, floor(min(y[, 2])))
    ymax <- max(w, ceiling(max(y[, 2])))

    out <- sapply(xmin:xmax, function(x, ...) {
        sapply(ymin:ymax, function(y, ...) {
            x <- ...[[1]]
            img <- ...[[2]]
            A <- ...[[3]]
            t <- ...[[4]]

            .x <- A %*% matrix(c(x,y), nrow=2) + t
            # .x <- round(.x)

            if (.x[1] < 1 | .x[1] > h | .x[2] < 1 | .x[2] > w) {
                # c(0, 0)
                0
            } else {
                # .x
                k <- interpolate(img, coords=.x, method="bilinear")
                k
                # img[.x[1], .x[2]]
            }

        }, c(list(x), ...))
    }, list(img, A, t))

    out <- as.Image(out)
    out
}


out <- sim_grid(c(500, 100), 10, FALSE)
t <- matrix(c(10, 0), nrow=2)
A <- matrix(c(1, 0.2, 0, 1), nrow=2, byrow=TRUE)
A <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)
display(out, "raster")
# out <- opt

tout <- sim_transform(out, A, t)
display(tout, "raster")

m <- matrix(c(1, -.5, 0, 0, 1, 0), nrow=3, ncol=2)
sim_transform(img, m[1:2, ], matrix(m[3, ], nrow=2))
affine(img, m)
display(.Last.value, "raster")
