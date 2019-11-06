

sample_coords <- function(dims, size=50) {
    # Width and heigh of image that is being sampled
    w <- dims[1]
    h <- dims[2]

    # Sampling and return as a matrix
    out <- matrix(
        c(sample.int(w, size=size, replace=TRUE),
            sample.int(n=h, size=size, replace=TRUE)),
        nrow = size, dimnames=list(NULL, c("x", "y")))

    out
}


intensities_from_coords <- function(img, coords, plot=FALSE, new_window=FALSE) {
    # check validity of img
    img <- .valid_img(img)

    # Dimensions of img
    dimx <- dim(img)[1]
    dimy <- dim(img)[2]

    if (!is.matrix(coords)) {
        stop("coords should be a matrix", call.=FALSE)}

    # Extract intensities (!!! dependence on EBImage)
    int <- apply(coords, 1, function(x) imageData(img)[x[1], x[2]])
    out <- cbind(coords, int)

    if (plot) {
        imageData(img) <- matrix(rep(0, dimx*dimy), nrow=dimx)
        int <- apply(out, 1, function(x) {
            imageData(img)[x[1], x[2]] <<- x[3]})

        if (new_window) quartz()
        display(img, method="raster")
    }

    out
}


coords <- sample_coords(dim(fixed), size=1000)
int <- intensities_from_coords(fixed, coords, plot=TRUE, new=TRUE)
