library("EBImage")
library("Cardinal")


.validate_overlay_images <- function(fixed, moving, alpha) {
    cls1 <- class(fixed)
    cls2 <- class(moving)
    if ( !cls1 == "Image" | !attr(cls1, "package") == "EBImage" |
         !cls2 == "Image" | !attr(cls2, "package") == "EBImage") {
        stop("Fixed and/or Moving images is not of 'Image' class
             of EBImage package", call.=FALSE)
    }

    if (alpha < 0 | alpha > 1) {
        stop("alpha has to be in [0,1] included", call.=FALSE)
    }
    fixed
}


overlay_images <- function(fixed, moving, alpha=0.25, new=FALSE) {
    .validate_overlay_images(fixed, moving, alpha)

    if (new){
        quartz()
    }
    display((alpha*fixed + (1-alpha)*moving), method="raster")
}



interpolate <- function(img, coords=c(h=1,w=1), method="bilinear") {
    h <- coords[[1]]
    w <- coords[[2]]

    h1 <- floor(h)
    h2 <- ceiling(h)
    w1 <- floor(w)
    w2 <- ceiling(w)

    int <- c()
    v <- matrix(c(h1,w1, h1,w2, h2,w1, h2,w2), ncol=2, byrow=TRUE)
    d <- apply(v, 1, function(x) {
        int <<- c(int, img[x[1], x[2]])
        dist(matrix(c(x, coords), ncol=2, byrow=TRUE), method="euclidean")
    })

    sum(int * d) / sum(d)
}

.k <- matrix(1:4, nrow=2, byrow=TRUE)
