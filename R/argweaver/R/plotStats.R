

##' Plot statistics trace vs iteration
##'
##' Given a data frame read in by readStatsFile, plot one of the
##' columns vs iteration number.
##' @param x A data.frame read in by readStatsFile. Alternatively, if
##' x is type character, can be the name of the stats file.
##' @param stat A character string naming the column of x to plot
##' @param add If FALSE, create a new plot. Otherwise add to current plot
##' @param xlim The x-axis range (only used when add==FALSE).
##' Taken from the range of data if not given.
##' @param ylim The y-axis range (only used when add==FALSE).
##' Taken from the range of data if not given
##' @param xlab The label for x-xais (only used when add=FALSE)
##' @param ylab The label for y-axis (only used when add=FALSE)
##' @param ... Other options to pass to plot function
##' @export
plotStats <- function(x, stat,
                     add=FALSE, xlim=NULL, ylim=NULL,
                     xlab="Iter", ylab=stat, ...) {
    
    if (length(stat) != 1 || class(stat) != "character")
        stop(stat, "should be a single character string")
    if (class(x) == "character") {
        if (length(x) != 1)
            stop("If x is a filename, should be length 1 (multiple files not supported)")
        x <- readStatsFile(x)
    }
    if (class(x) != "data.frame")
        stop("x should be data.frame")
    if (!is.element(stat, names(x)))
        stop("Did not find column ", stat, " in x")
    if (!is.element("iter", names(x)))
        stop("Did not find \"iter\" column in x")
    if (!add) {
        if (!is.null(xlim)) {
            xlim <- as.numeric(xlim)
            if (length(xlim) != 2)
                stop("Invalid xlim")
        } else {
            xlim <- range(x$iter)
        }
        if (!is.null(ylim)) {
            ylim <- as.numeric(ylim)
            if (length(ylim) != 2)
                stop("Invalid ylim")
        } else {
            ylim <- range(x[,stat])
        }
        plot(x$iter, x[,stat], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
             ...)
    } else {
        points(x$iter, x[,stat], ...)
    }
    invisible(NULL)
}
