##' Overlay several histograms
##'
##' x should be a list of data frames with the same structure.
##' Will plot a histogram of the column named by "stat"
##' for each item in x
##' @param x list of data frames. Each should have a column named stat
##' @param stat gives the column to plot. Should be single character string.
##' @param col A list or character vectors of colors to use for each histogram.
##' If NULL, will use a rainbow color, but each color will be given transparency
##' value alpha=0.5
##' @param xlab Label for x axis
##' @param breaks Passed on to hist function. Can be a number of breaks, or a numeric
##' vector of break intervals.
##' @param main Title for plot
##' @param xlim limits for x-axis
##' @param ... Passed to "hist" function
##' @export
multiHist <- function(x, stat, col=NULL, xlab=stat, breaks=50, main="", xlim=NULL, ...) {
    if (class(x) != "list")
        stop("x should be a list")
    defaultColors <- (is.null(col))
    for (i in 1:length(x)) {
        if (!is.element(stat, names(x[[i]])))
            stop("Did not find stat ", stat, " in x[[",i,"]]")
        if (i==1) {
            alldata <- x[[i]][,stat]
        } else alldata <- c(alldata, x[[i]][,stat])
    }
    tmp <- hist(alldata, plot=FALSE, breaks=breaks, ...)
    breaks <- tmp$breaks
    for (i in 1:length(x)) {
        tmp <- hist(x[[i]][,stat], breaks=breaks, plot=FALSE, ...)
        if (i==1) ylim <- range(tmp$counts) else ylim <- range(c(ylim, tmp$counts))
    }
    for (i in 1:length(x)) {
        if (defaultColors) {
            useCol <- rainbow(length(x), alpha=0.5)[i]
        } else if (is.list(col)) {
            useCol <- col[[names(x)[i]]]
        } else useCol <- col[i]
        if (is.null(xlim)) xlim=range(alldata)
        hist(x[[i]][,stat], breaks=breaks, ylim=ylim, xlim=xlim, col=useCol,
             add=(i > 1), xlab=xlab, main=main, ...)
    }
    invisible(NULL)
}
