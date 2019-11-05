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
