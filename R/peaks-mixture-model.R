library(mclust)
library(ggplot2)


# Split well into groups by bead ID, then apply 'FUN' to each group in turn.
applyWell <- function(well, FUN=identity)
{
    tapply(well[ ,'RP1'], well[ ,'RID'], FUN)
}

# Apply 'FUN' to a particular bead in all wells, one well at a time.
applyBead <- function(wells, bead, FUN=identity)
{
    lapply(wells, function(well) FUN(well[well[ ,'RID'] == bead, 'RP1']))
}

# Print summary of peaks and median for the given well.  Call
#
#   lapply(wells, peakSummary)
#
# to print a summary for a list of wells.
peakSummary <- function(well)
{
    data.frame(best=applyWell(well,FUN=bestPeak),
               left=applyWell(well,FUN=leftPeak),
               median=applyWell(well,FUN=median),
               right=applyWell(well,FUN=rightPeak))
}

# Plot each bead density in well and mark left&right peaks as well as the
# median.
plotPeaks <- function(well)
{
    devAskNewPage(TRUE)
    applyWell(well, plotPeak)
    devAskNewPage(FALSE)
}

# Plot density for 'xlin' and mark the left&right peaks as well as the median.
plotPeak <- function(xlin)
{
    p  <- suppressWarnings(estimatePeaks(xlin))
    df <- as.data.frame(p)
    print(df)
    if (any(is.na(df)))
        stop('Peak estimation failed, skipping...')

    # Detect which peak to pick:
    # Take the one with smallest variability, unless it has a proportion
    # smaller than 'eps'.
    eps <- 0.25
    l <- p$var[1] < p$var[2]
    if (l && p$prop[1] < eps)  l <- FALSE
    if (!l && p$prop[2] < eps) l <- TRUE

    hist <- qplot(xlin, geom="density", log="x")
    print(hist +
          geom_histogram(aes(y = ..density..), alpha = 0.1) +
          geom_vline(aes(xintercept=p$mean[1]), col="red",
                     linetype=ifelse(l, "solid", "dashed")) +
          geom_vline(aes(xintercept=p$mean[2]), col="blue",
                     linetype=ifelse(l, "dashed", "solid")) +
          geom_vline(aes(xintercept=median(xlin)), col="green",
                     linetype="dashed"))
}

# Use two-Gaussian mixture model to detect the peaks (i.e. means of the two
# Gaussians) of 'xlin'.  The data 'xlin' is log-transformed before applying the
# mixture model.  Hence for the estimate to be somewhat accurate the data
# should have two near-Gaussian bumps after being log-transformed.
#
# The return value is a list which describes the proportion, mean, and variance
# of the fitted Gaussians.
estimatePeaks <- function(xlin)
{
    # The data is assumed to be linear.  Take logarithm first so that the data
    # is closer to Gaussian.
    x <- log(xlin[!is.na(xlin) & xlin>0])

    # Need to supply initial guess as to which group each sample belongs to.
    # Our guess is:
    #   - Group 1, if sample <= median
    #   - Group 2, if sample >  median
    guess <- unmap(x > median(x))

    # Perform EM algorithm to estimate parameters
    #   "E" - equal variance
    #   "V" - variable variance
    ms <- me("V", x, guess)
    es <- em("V", x, ms$parameters)

    # Mean is converted back to linear.
    par <- es$parameters
    list(prop=par$pro, mean=exp(par$mean), var=par$variance$sigmasq)
}

# Pick out left peak from two-Gaussian mixture model.
leftPeak <- function(xlin)
{
    p <- suppressWarnings(estimatePeaks(xlin))
    p$mean[1]
}

# Pick out right peak from two-Gaussian mixture model.
rightPeak <- function(xlin)
{
    p <- suppressWarnings(estimatePeaks(xlin))
    p$mean[2]
}

# Pick out "best" peak from two-Gaussian mixture model.
bestPeak <- function(xlin, eps=0.25)
{
    p  <- suppressWarnings(estimatePeaks(xlin))
    df <- as.data.frame(p)
    if (any(is.na(df)))
        return(NA)

    # Detect which peak to pick:
    # Take the one with smallest variability, unless it has a proportion
    # smaller than 'eps'.
    l <- p$var[1] < p$var[2]
    if (l && p$prop[1] < eps)  l <- FALSE
    if (!l && p$prop[2] < eps) l <- TRUE

    if (l) p$mean[1] else p$mean[2]
}
