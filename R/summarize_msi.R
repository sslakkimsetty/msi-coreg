library("Cardinal")
library("EBImage")

dpath1 <- "/Users/sai/Documents/00-NEU/2-Ind-Study/data"
dpath2 <- "S042"
dpath3 <- ""
dpath <- file.path(dpath1, dpath2, dpath3)

msifname <- "S042-processed.imzML"

# Read MSI data
msi2 <- readMSIData(file.path(dpath, msifname))
plot(msi2, pixel=1)

msi <- process(normalize(msi, method="tic"))
msi <- process(smoothSignal(msi, method="gaussian"))
msi <- process(reduceBaseline(msi, method="median", blocks=9))
msi <- process(peakPick(msi, method="adaptive", SNR=3))

msi2 <- preprocess_moving(msi2)
all(iData(msi) >= 0)



# HR2MSI DATASET
dpath1 <- "/Users/sai/Documents/00-NEU/2-Ind-Study/data"
dpath2 <- "HR2MSI-mouse-unary-bladder"
dpath3 <- ""
dpath <- file.path(dpath1, dpath2, dpath3)

msifname <- "HR2MSI mouse urinary bladder S096.imzML"
hr2msi <- readMSIData(file.path(dpath, msifname))

optfame <- "HR2MSI mouse urinary bladder S096 - optical image.tiff"
hr2opt <- readImage(file.path(dpath, optfname))

opt <- preprocess_fixed(hr2opt)
msi <- preprocess_moving(hr2msi)

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

# Run spatial DGMM with top m/z values (DGMM fails!!!)
a <- features(msi, mz=c(559.06, 633))
dgmm <- spatialDGMM(msi[tissue][a, ], r=1, k=3, method="adaptive", init="kmeans",
    iter.max=1000, tol=1e-11)



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

a <- out[[2]]

a <- as.vector(t(a))

image(msi, mz=a[1:40], superpose=TRUE, key=FALSE)

image(out[[1]][[1]], values="probability")

b <- out[[1]][[1]]
str(b)

head(as.integer(b@resultData@listData[[1]]$class))
str(b@resultData@listData[[1]]$class)

b@resultData@listData[[1]]$class <- b@resultData@listData[[1]]$class + 1

lev <- levels(b@resultData@listData[[1]]$class)

levels(b@resultData@listData[[1]]$class) <- c(lev, "d")

b@resultData@listData[[1]]$class[[1]] <- "d"
