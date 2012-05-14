# Plot binned data.
#
# The parameter 'x' is the bin labels (i.e. categories, factors) and 'y' is a
# vector of the same length as 'x'.  Both variable names may be columns in the
# data frame 'data' (if specified).  A line connecting the median of each bin
# is drawn in color 'medianColor'.  The interquartile range of the bins are
# shaded in color 'quartileColor' (don't shade if this is NULL).  Finally, the
# range between the min and max of each bin is shaded in color 'rangeColor'
# (don't shade if this is NULL).
#
# Examples:
#
# plotBins(rep(1:10,10), rnorm(100))
#
# Example on how to plot male plasma samples for antibody 1901 and sort them
# into 5 year bins (this assumes df$age, df$gender, df$type, and df$Ab_1891
# exist):
#
#   df <- loadAGE()
#   plotBins(bins(age,width=5), Ab_1891,
#            data=subset(df, gender == 'male' & type == 'plasma'),
#            main="A protein profile", xlab="Age" ylim=c(185,235))
#
plotBins <- function(x, y,
                     data = data.frame(),
                     medianColor = 'blue',
                     quartileColor = rgb(0.5,0.5,0.5,0.3),
                     rangeColor = rgb(0.5,0.5,0.5,0.1),
                     xlab = deparse(substitute(x)),
                     ylab = deparse(substitute(y)), ...)
{
    # Evaluate the 'x' and 'y' arguments inside the database 'data' so that
    # they can refer to variables inside 'data'.
    attach(data)
    X <- x
    Y <- y
    detach()

    # Scatter plot data points ('plot' will do box-and-whisker plot if 'X' is a
    # factor, hence the use of 'plot.default').
    plot.default(X, Y, xlab=xlab, ylab=ylab, ...)

    # Determine quantiles
    qs <- tapply(Y, X, quantile, na.rm=TRUE)
    # Figure out which x coordinates to use
    qx <- if (is.factor(X)) 1:length(qs) else names(qs)
    # Convert quantiles into matrix whose rows are: min,q1,median,q3,max
    qm <- matrix(unlist(qs), nrow=5)

    # Draw shaded area between min and max
    if (!is.null(rangeColor))
        polygon(c(qx,rev(qx)), c(qm[1,],rev(qm[5,])), col=rangeColor)

    # Draw shaded area between 1st and 3rd quartiles
    if (!is.null(quartileColor))
        polygon(c(qx,rev(qx)), c(qm[2,],rev(qm[4,])), col=quartileColor)

    # Draw line on median
    lines(qx, qm[3,], col=medianColor, lw=3)
}

bins <- function(x, n=5)
{
    bs <- seq(from=n/2, to=max(x)+n, by=n)
    ls <- seq(from=n, by=n, length.out=length(bs)-1)
    a  <- cut(x, breaks=bs, labels=ls)
    as.numeric(levels(a))[a]
}

randomProfile <- function(binCount, binSize)
{
    data.frame(x=rep(1:binCount, binSize), y=rnorm(binCount*binSize))
}

# Compute signature scores.  The parameter 'cats' contains the bin labels, and
# 'pts' is a vector of the same length as 'cats'.
#
# NOTE! This function requires the 'sigscore' tool to be in the $PATH.
#
sigscore <- function(cats, pts)
{
    if (is.vector(pts) && length(cats) == length(pts))
        pts <- as.matrix(pts)

    if (length(cats) != dim(pts)[1])
        stop('Number of categories does not match rows of data')

    m <- rbind(cats, t(pts))

    message('Writing temporary CSV file...')
    fname <- tempfile()
    csv <- write.table(m, file=fname, row.names=FALSE, col.names=FALSE,
                       na="", sep=",")

    message('Computing signature scores...')
    z <- system(paste('sigscore',fname), intern=TRUE)
    system(paste('rm',fname))
    as.numeric(z)
}
