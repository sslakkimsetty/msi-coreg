library("Cardinal")

dpath1 <- "/Users/sai/Documents/00-NEU/2-Ind-Study/data"
dpath2 <- "S042"
dpath3 <- ""
dpath <- file.path(dpath1, dpath2, dpath3)

msifname <- "S042-continuous.imzML"

# Read MSI data
msi <- readMSIData(file.path(dpath, msifname))


# Simulate MSI data
register(SerialParam())
set.seed(2020)

msi <- simulateImage(preset=1, npeaks=50, nruns=1, baseline=2)


preprocess_moving <- function(moving,
    process=c("normalize", "smooth", "baseline", "pick", "align", "filter"),
    seed=210) {

    if ("normalize" %in% process) {
        print("Normalizing...")
        moving <- normalize(moving, method="tic")
    }

    if ("smooth" %in% process) {
        print("Smoothing...")
        moving <- smoothSignal(moving, method="gaussian", window=9)
    }

    if ("baseline" %in% process) {
        print("Correcting baseline...")
        moving <- reduceBaseline(moving, method="median", blocks=50)
    }

    if ("pick" %in% process) {
        print("Peak picking...")
        moving <- peakPick(moving, method="adaptive", SNR=3)
    }

    if ("align" %in% process) {
        if ( !("pick" %in% process) ) {
            stop("Peak picking should also be included!", call.=FALSE)
        }
        print("Peak aligning...")
        moving <- peakAlign(moving, method="diff")
    }

    if ("filter" %in% process) {
        if ( !all(c("pick", "align") %in% process) ) {
            stop("Peak pick or align missing and should be included!",
                call.=FALSE)
        }
        print("Peak filtering...")
        # moving <- peakFilter(moving, freq.min=0.2)
    }

    moving <- process(moving)

    moving
}

msi2 <- preprocess_moving(msi)
plot(msi2, pixel=2)

