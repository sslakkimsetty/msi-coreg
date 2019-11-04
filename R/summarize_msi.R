library("Cardinal")
library("EBImage")

# HR2MSI dataset
dpath1 <- "/Users/sai/Documents/00-NEU/2-Ind-Study/data"
dpath2 <- "HR2MSI-mouse-unary-bladder"
dpath3 <- ""
dpath <- file.path(dpath1, dpath2, dpath3)

msifname <- "HR2MSI mouse urinary bladder S096.imzML"
hr2msi <- readMSIData(file.path(dpath, msifname))

optfname <- "HR2MSI mouse urinary bladder S096 - optical image.tiff"
hr2opt <- readImage(file.path(dpath, optfname))

opt <- preprocess_fixed(hr2opt) # processed optical image
msi <- preprocess_moving(hr2msi) # processed msi image

# Run spatial shrunken centroids and spatial DGMM to find optimal m/z
RNGkind("L'Ecuyer-CMRG")
set.seed(2019)
ssc <- spatialShrunkenCentroids(hr2msi, r=1, k=3,
    s=c(0,3), method="gaussian", BPPARAM=MulticoreParam())

ssc2 <- spatialShrunkenCentroids(hr2msi, r=1, k=3,
    s=c(6,9), method="gaussian", BPPARAM=MulticoreParam())

ssc_sim <- spatialShrunkenCentroids(msi, r=1, k=3,
    s=c(6,9), method="gaussian", BPPARAM=MulticoreParam())

image(ssc2, model=list(s=9), values="probability")
image(ssc_sim, model=list(s=9), values="probability")

top_ssc2 <- topFeatures(ssc2, model=list(s=6))
top_ssc_sim <- topFeatures(ssc_sim, model=list(s=6))

image(hr2msi, mz=633)
image(msi, mz=633)

ssc3 <- spatialShrunkenCentroids(msi[tissue], r=1, k=3,
    s=c(6,9), method="gaussian", BPPARAM=MulticoreParam())
image(ssc3, model=list(s=9), values="probability", key=FALSE)


.validate_summarize_input <- function(moving, r, k, s) {
    moving <- .validate_moving_input(moving)

    if ( !all(c(typeof(r), typeof(k), typeof(s)) == "double") |
        !all(c(class(r), class(k), class(s)) == "numeric") ) {
        stop("r, k, and s must be of numeric class and double type",
            call. = FALSE)
    }

    if ( any(k >= 6) ) {
        warning("k above 5 is not recommended")
    }

    if ( any(r > 3) ) {
        warning("r above 3 is not recommended")
    }

    moving
}


.compute_ssc <- function(moving, r, k, s) {
    ssc <- spatialShrunkenCentroids(moving, r=r, k=k, s=s,
        method="adaptive", BPPARAM=MulticoreParam())

    top <- topFeatures(ssc)$mz

    return(list(ssc, top))
}


summarize_moving_ssc <- function(moving, r=c(1,2),
    k=c(3,4,5), s=c(0,3,6,9), seed=210) {

    RKS <- sapply(list(r, k, s), function(x) length(x))
    size <- cumprod(RKS)[3]

    RNGkind("L'Ecuyer-CMRG")
    moving <- .validate_summarize_input(moving, r, k, s)

    top <- matrix(rep(0, 10*size), nrow=10)
    colnames(top) <- rep("name", size)

    ssc_list <- c()

    count <- 1

    for (.s in s) {
        for (.k in k) {
            for (.r in r) {
                out <- .compute_ssc(moving, .r, .k, .s)
                ssc_list <- c(ssc_list, out[[1]])
                top[, count] <- out[[2]]

                colnames(top)[count] <- paste0("s=", .s, ", k=", .k, ", r=", .r)
                count <- count + 1
            }
        }
    }

    return(list(ssc_list, top))

}


out <- summarize_moving_ssc(msi, r=c(1), k=c(3,4), s=c(0,3,6,9))


# Build composite from ssc images
fac <- b@resultData@listData[[1]]$class
dims(msi)

fac <- as.numeric(levels(fac))[fac]
fac <- matrix(fac, nrow=260)

fac_img <- Image(fac)
fac_img <- preprocess_fixed(fac_img)
display(fac_img, "raster")


fac_l <- lapply(out[[1]], function(x) x@resultData@listData[[1]]$class)
fac_l <- lapply(fac_l, function(x) as.numeric(levels(x))[x])

build_image <- function(fac, nrow) {
    fac <- matrix(fac, nrow=nrow)
    fac_img <- Image(fac)
    fac_img <- preprocess_fixed(fac_img)
    fac_img
}

fac_img_l <- lapply(fac_l, function(x) build_image(x, nrow=260))

sum <- fac_img_l[[1]]

for (i in c(2, 3, 4, 6, 7, 8)) {
    sum <- sum + fac_img_l[[i]]
}

fac_avg <- sum / length(fac_img_l)
display(fac_avg, "raster")


# Build composite from top mz values
top_mz = as.vector(out[[2]])
image(msi, mz=top_mz, superpose=TRUE, key=FALSE)

# Build composite image from multiple mz values
composite_img_mz <- function(msi, top_mz) {
    N <- length(top_mz)

    img_l <- lapply(top_mz,
        function(x) as.vector(image(msi, mz=x)[[1]][[1]][[1]][[3]]))
    img_l <- lapply(img_l, function(x) build_image(x, nrow=260))

    composite <- img_l[[1]]
    for (i in 2:N) {
        composite <- composite + img_l[[i]]
    }

    composite <- preprocess_fixed(composite)
}

x <- composite_img_mz(msi, top_mz)
.x <- thresh(x, w=5, h=5, offset=0.0001)
str(x)
display(.x, "raster")

.fac_avg <- thresh(fac_avg, w=1, h=1, offset=0.02)
display((.fac_avg + 3*fac_avg)/4, "raster")
