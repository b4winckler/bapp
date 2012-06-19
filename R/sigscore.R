# Compute signature scores.  The parameter 'cats' is a vector of bin labels,
# and 'pts' is a matrix with as many rows as there are elements in 'cats'.
#
# Computing the signature score can take a long time and for this reason a
# progress bar is displayed during the computation.  Set 'quiet=TRUE' to
# suppress the progress bar.  Use fewer categories to speed up the computation
# (the 'bins()' function can be used for this).
#
# NOTE! This function requires the 'sigscore' tool to be in the $PATH.
#
sigscore <- function(cats, pts, quiet=FALSE)
{
    if (is.vector(pts))
        pts <- as.matrix(pts)

    d <- dim(pts)
    if (length(cats) != d[1])
        stop('Number of categories does not match rows of data')

    m <- rbind(cats, t(pts))

    if (!quiet)
        message(paste('Computing', d[2], 'signature scores...'))
    fname <- tempfile()
    csv <- write.table(m, file=fname, row.names=FALSE, col.names=FALSE,
                       na="", sep=",")

    z <- system(paste('sigscore',fname), intern=TRUE, ignore.stderr=quiet)
    system(paste('rm',fname))

    vals <- as.numeric(z)
    if (!is.null(colnames(pts)))
        names(vals) <- colnames(pts)

    vals
}

