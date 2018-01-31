##' plot bigWig (a.k.a. bedGraph-style) data
##'
##' @param x a data frame with at least three columns with names
##' chrom, chromStart, chromEnd
##' @param y a numeric vector whose length=nrow(x) giving one score per
##' region represented in x
##' @param col The color of the line to draw
##' @param alpha Transparency factor (0 = transparent, 1 = solid)
##' @param y0 A lower limit to y (should be numeric vector, same length as x and y)
##' @param y1 An upper limit to y (same format as y0; must be provided if y0 is provided)
##' @param add If TRUE, add to current plot
##' @param xlim The x-axis limits
##' @param chrom If given plot only elements of x for each x$chrom==chrom
##' @param ylab label for y-axis
##' @param xlab label for x-axis
##' @param ylim limits for y-axis. If not given will use range(y)
##' @param ... Other arguments to be passed to plot functions
##'
##' If y0 and y1 are given, then will shade in regions between y0 and y1 in the same
##' color as the line drawn, but made transparent by factor "alpha'
##' #TODO: have it take an optional low/high value, as well as optional alpha
##' @export
plotBw <- function(x, y, y0=NULL, y1=NULL, col="black", alpha=0.2,
                   add=FALSE, xlim=NULL, chrom=NULL, ylab="Generations",
                   xlab="Coordinate", ylim=NULL, ...) {
    if (!is.null(y0) || !is.null(y1)) {
        if (is.null(y0) || is.null(y1))
            stop("plotBw: need to supply both y0 and y1, or neither")
        if (length(y0) != length(y) || length(y1) != length(y))
            stop("plotBw: y0 and y1 length should match y length")
    }
    x$chrom <- as.character(x$chrom)
    if (!is.null(chrom)) {
        f <- is.element(x$chrom, chrom)
        x <- x[f,]
        y <- y[f]
    }
    allChroms <- unique(x$chrom)
    numChrom <- length(allChroms)
    chromLen <- numeric(numChrom)

    for (i in 1:numChrom)
        chromLen[i] <- max(x[x$chrom == allChroms[i],"chromEnd"])
    if (numChrom > 1) {
        for (i in 2:numChrom) {
            f <- x$chrom == allChroms[i]
            addnum <- sum(chromLen[1:(i-1)])
            x[f,"chromStart"] <- x[f,"chromStart"] + addnum
            x[f,"chromEnd"] <- x[f,"chromEnd"] + addnum
        }
    }

    if (!add) {
        if (is.null(xlim)) xlim <- range(c(x$chromStart, x$chromEnd))
        if (is.null(ylim)) ylim <- range(y)
        plot(c(0), c(0), xlab=xlab, ylab=ylab, type="n", xlim=xlim, ylim=ylim, xaxt="n", col=col, ...)
    }
    midChrom <- numeric(numChrom)
    for (i in 1:numChrom) {
        f <- x$chrom == allChroms[i]
        midChrom[i] <- mean(c(min(x[f,"chromStart"]), max(x[f,"chromEnd"])))
        if (!is.null(y0))
            rect(xleft=x$chromStart[f], xright=x$chromEnd[f],
                 ybottom=y0[f], ytop=y1[f], col=makeTransparent(col, alpha),
                 border=NA)
        segments(x0=x$chromStart[f], x1=x$chromEnd[f], y0=y[f], col=col, ...)
        segments(x0=x[f,]$chromStart[-1], y0=y[f][-1], y1=y[f][-nrow(x)], col=col, ...)
    }
    if (!add) {
        if (numChrom >= 2) {
            for (i in 2:numChrom) {
                abline(v=sum(chromLen[1:(i-1)]), lty=3)
            }
            mtext(allChroms, side=1, line=0.5, at=midChrom, cex=0.7, las=2)
        } else axis(1)
    }
    invisible(NULL)
}
