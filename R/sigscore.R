# Compute signature scores.  The parameter 'cats' contains the bin labels, and
# 'pts' is a vector of the same length as 'cats'.
#
# NOTE! This function requires the 'sigscore' tool to be in the $PATH.
#
sigscore <- function(cats, pts)
{
    if (is.vector(pts))
        pts <- as.matrix(pts)

    if (length(cats) != dim(pts)[1])
        stop('Number of categories does not match rows of data')

    m <- rbind(cats, t(pts))

    message('Computing signature scores...')
    fname <- tempfile()
    csv <- write.table(m, file=fname, row.names=FALSE, col.names=FALSE,
                       na="", sep=",")

    z <- system(paste('sigscore',fname), intern=TRUE)
    system(paste('rm',fname))
    as.numeric(z)
}

