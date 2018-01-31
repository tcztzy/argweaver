##' Read in a LeafTrace file
##' @param file The leaf trace file (produced by the script arg-layout)
##' @return A data frame that can be set to plotLeafTrace
##' @export
readLeafTrace <- function(file) {
    x <- read.table(file, stringsAsFactors=FALSE, header=FALSE)
    names(x)[1:3] <- c("chrom", "chromStart", "chromEnd")
    x
}

##' Create a leafTrace plot
##'
##' The leaftrace draws a horizontal line for every leaf across the ARG,
##' such that the distance between adjacent lines reflects the distance between
##' the leaves corresponding to those lines. There is an randomness to the
##' leaf order and the distance between non-adjacent lines cannot be ascertained
##' 
##' @param file A leaf trace file (produced by the script arg-layout)
##' @param x (alternative to file) A leaf trace data frame read in by readLeafTrace
##' @param col Either a character string giving the color of the plot, or a
##' list indexed by leaf or diploid individual names giving the color for
##' each leaf/individual. Will assume that leafs are named with the
##' convention <ind>_1 and <ind>_2 for the two diploid lineages
##' @param xlim The range for x coordinates
##' @param ylim The range for y coordinates
##' @param xlab Label for the x-axis
##' @param ylab Label for the y-axis
##' @param lwd Line width
##' @param add If TRUE, add plot to current plot
##' @param subset If NULL, plot all leafs. Otherwise only plot individuals/leafs
##' in this character vector
##' @param ... Additional objects passed to the plot function
##' @export
plotLeafTrace <- function(file=NULL, x=NULL, col="black", xlim=NULL, ylim=NULL, xlab="Coordinate", ylab="", add=FALSE, subset=NULL, lwd=1, ...) {
    if (is.null(file) && is.null(x)) {
        stop("Need to provide either a layout file (file) or a layout table (x)")
    }
    if (!is.null(file))
        x <- readLeafTrace(file)
    if (is.null(xlim))
        xlim <- range(c(x$chromStart, x$chromEnd))
    if (!is.null(subset))
        subset <- c(subset, paste0(subset, "_1"), paste0(subset, "_2"))
    if (is.list(col)) {
        currNames <- names(col)
        col[sprintf("%s_1", currNames)] <- col[currNames]
        col[sprintf("%s_2", currNames)] <- col[currNames]
    }
    nameCols <- seq(from=4, to=ncol(x), by=2)
    yCols <- nameCols+1
    if (is.null(ylim))
        ylim <- range(x[,yCols])
    if (!add)
        plot(c(0),c(0), type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    for (i in nameCols) {
        ind <- unique(x[,i])
        if (length(ind) != 1) stop("got different inds in same col")
        if (!is.null(subset) && !is.element(ind, subset)) next
        if (is.list(col)) {
            color <- as.character(col[x[,i]])
            color2 <- color[-1]
        } else {
            color <- col
            color2 <- col
        }
        segments(x0=x$chromStart, x1=x$chromEnd, y0=x[,i+1], col=color, lwd=lwd)
        segments(x0=x$chromStart[-1], y0=x[,i+1][-1], y1=x[,i+1][-nrow(x)], col=color2, lwd=lwd)
    }
    invisible(NULL)
}
