## ----setup, include = FALSE------------------------------------------------
getCaps <- function ()
{
    label <- knitr::opts_current$get("label")
    txts <- list(
        figcover=c(
            "Variation in the sharpness of the peaks using \\texttt{trim}",
            "attribute."
        ),
        figmnase=c(
            "Toy example of MNase biase correction. Random nucleosomal and",
            "control reads have been generated using \\texttt{synteticNucMap}",
            "function and corrected using \\texttt{controlCorrect}."
        ),
        fignoise=c(
            "Original intensities from tiling array experiment. Smoothing",
            "using a sliding window of variable length (0, 20, 50 and 100 bp)",
            "is presented."
        ),
        figfft=c(
            "Power spectrum of the example Tiling Array data, percentile 1",
            "marked with a dashed line."
        ),
        figfft2=c(
            "Filtering in Tiling Array (up, blue) (1\\% comp.) and NGS (down",
            "red) (2\\% comp.)."
        ),
        figpeak=c(
            "Output of \\texttt{plotPeaks} function. Peaks are spotted in red",
            "and detection threshold marked with an horitzontal line."
        ),
        figpeak2="\\texttt{plotPeaks} function with \\texttt{score=TRUE}.",
        figpeak3=c(
            "\\texttt{plotPeaks} output with \\texttt{score=TRUE} and",
            "\\texttt{width=140}."
        ),
        figranges=c(
            "Simple example of ranges manipulation to plot fuzzy nucleosomes."
        ),
        figsyn=c(
            "Example synthetic coverage map of 90 well-positioned (100-10)",
            "and 20 fuzzy nucleosomes."
        )
    )
    caps <- lapply(txts, paste, collapse=" ")
    return (caps[[label]])
}

knitr::opts_chunk$set(
    collapse   = TRUE,
    comment    = "#>",
    message    = FALSE,
    warning    = FALSE,
    fig.width  = 6,
    fig.height = 2,
    fig.wide   = TRUE
)


## ----processtilling, eval=FALSE--------------------------------------------
#  processTilingArray(data, exprName, chrPattern, inferLen=50)

## ----loadTA----------------------------------------------------------------
library(nucleR)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
data(nucleosome_tiling)
head(nucleosome_tiling, n=25)

## ----readBAM---------------------------------------------------------------
sample.file <- system.file("extdata", "cellCycleM_chrII_5000-25000.bam",
    package="nucleR")
reads <- readBAM(sample.file, type="paired")
head(reads)

## ----import----------------------------------------------------------------
data(nucleosome_htseq)
class(nucleosome_htseq)
nucleosome_htseq

## ----processReads----------------------------------------------------------
# Process the paired end reads, but discard those with length > 200
reads_pair <- processReads(nucleosome_htseq, type="paired", fragmentLen=200)

# Process the reads, but now trim each read to 40bp around the dyad
reads_trim <- processReads(nucleosome_htseq, type="paired", fragmentLen=200,
    trim=40)

## ----coverage--------------------------------------------------------------
# Calculate the coverage, directly in reads per million (r.p.m)
cover_pair <- coverage.rpm(reads_pair)
cover_trim <- coverage.rpm(reads_trim)

## ----figcover, echo=FALSE, fig.cap=getCaps(), fig.width=5------------------
# Compare both coverages
t1 <- as.vector(cover_pair[[1]])[1:2000]
t2 <- as.vector(cover_trim[[1]])[1:2000]
t1 <- (t1 - min(t1)) / max(t1 - min(t1)) # Normalization
t2 <- (t2 - min(t2)) / max(t2 - min(t2)) # Normalization
plot_data <- rbind(
    data.frame(
        x=seq_along(t1),
        y=t1,
        coverage="original"
    ),
    data.frame(
        x=seq_along(t1),
        y=t2,
        coverage="trimmed"
    )
)
ggplot(plot_data, aes(x=x, y=y)) +
    geom_line(aes(color=coverage)) +
    xlab("position") +
    ylab("norm coverage")

## ----figmnase, echo=c(1:5), fig.cap=getCaps(), fig.width=5-----------------
# Toy example
map <- syntheticNucMap(as.ratio=TRUE, wp.num=50, fuz.num=25)
exp <- coverage(map$syn.reads)
ctr <- coverage(map$ctr.reads)
corrected <- controlCorrection(exp, ctr)

plot_data <- rbind(
    data.frame(
        x=seq_along(exp),
        y=as.vector(exp),
        coverage="normal"
    ),
    data.frame(
        x=seq_along(corrected),
        y=as.vector(corrected),
        coverage="corrected"
    )
)
ggplot(plot_data, aes(x=x, y=y)) +
    geom_line(aes(color=coverage)) +
    xlab("position") +
    ylab("coverage")

## ----fignoise, echo=FALSE, fig.cap=getCaps(), fig.height=1.5---------------
windowFilter <- function (x, w) {
    if (missing(w)) {
        return(x)
    } else {
        y <- filter(x, rep(1, w)/w)
        return(as.vector(y))
    }
}

mkEntry <- function (x, i, w, lab) {
    if (missing(lab)) {
        lab <- sprintf("slinding w. %i bp", w)
    }
    df <- data.frame(x=i, y=windowFilter(x[i], w), lab=lab)
    df[!is.na(df[, "y"]), ]
}

i <- 1:2000
plot_data <- rbind(
    mkEntry(nucleosome_tiling, i, 1, "original"),
    mkEntry(nucleosome_tiling, i, 20),
    mkEntry(nucleosome_tiling, i, 50),
    mkEntry(nucleosome_tiling, i, 100)
)

ggplot(plot_data, aes(x=x, y=y)) +
    geom_line() +
    facet_grid(.~lab) +
    xlab("position") +
    ylab("intensity")

## ----figfft, fig.cap=getCaps(), fig.width=5--------------------------------
fft_ta <- filterFFT(nucleosome_tiling, pcKeepComp=0.01, showPowerSpec=TRUE)

## ----figfft2, echo=FALSE, fig.cap=getCaps(), fig.width=5, fig.height=4-----
i <- 1:3000
tiling_raw <- nucleosome_tiling[i]
tiling_fft <- filterFFT(tiling_raw, pcKeepComp=0.01)
htseq_raw <- as.vector(cover_trim[[1]])[i]
htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.02)

plot_data <- rbind(
    data.frame(x=i, y=tiling_raw, lab="intensity", treatment="raw"),
    data.frame(x=i, y=tiling_fft, lab="intensity", treatment="filtered"),
    data.frame(x=i, y=htseq_raw,  lab="coverage",  treatment="raw"),
    data.frame(x=i, y=htseq_fft,  lab="coverage",  treatment="filtered")
)

ggplot(plot_data, aes(x=x, y=y)) +
    geom_line(aes(color=treatment)) +
    facet_grid(lab~., scales="free_y") +
    ylab("") +
    xlab("position")

## ----corfft----------------------------------------------------------------
tiling_raw <- nucleosome_tiling
tiling_fft <- filterFFT(tiling_raw, pcKeepComp=0.01)
htseq_raw <- as.vector(cover_trim[[1]])
htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.02)

cor(tiling_raw, tiling_fft, use="complete.obs")
cor(htseq_raw, htseq_fft, use="complete.obs")

## ----figpeak, fig.cap=getCaps()--------------------------------------------
peaks <- peakDetection(htseq_fft, threshold="25%", score=FALSE)
peaks

plotPeaks(peaks, htseq_fft, threshold="25%", ylab="coverage")

## ----figpeak2, fig.cap=getCaps()-------------------------------------------
peaks <- peakDetection(htseq_fft, threshold="25%", score=TRUE)
head(peaks)

plotPeaks(peaks, htseq_fft, threshold="25%")

## ----figpeak3, fig.cap=getCaps()-------------------------------------------
peaks <- peakDetection(htseq_fft, threshold="25%", score=TRUE, width=140)
peaks

plotPeaks(peaks, htseq_fft, threshold="25%")

## ----figranges, echo=-6, fig.cap=getCaps()---------------------------------
nuc_calls <- ranges(peaks[peaks$score > 0.1, ])
red_calls <- reduce(nuc_calls)
red_class <- RangedData(red_calls, isFuzzy=width(red_calls) > 140)
red_class

plotPeaks(red_calls, htseq_fft, threshold="25%")

## ----figsyn, fig.cap=getCaps(), fig.height=3-------------------------------
syn <- syntheticNucMap(wp.num=100, wp.del=10, wp.var=30, fuz.num=20,
    fuz.var=50, max.cover=20, nuc.len=147, lin.len=20, rnd.seed=1,
    as.ratio=TRUE, show.plot=TRUE)

