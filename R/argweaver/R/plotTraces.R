
##' Plot ARGweaver statistics vs iteration
##'
##' Plots the statistics in stat files output by arg-sample vs iteration.
##' By default plots all statistics given in the stat file (as well
##' as any extra created by readStatFiles).
##' @param files A character vector of file names giving the stat files to plot
##' @param dir (alternative to files; ignored if files is given). Plot all stat files found in this directory.
##' @param outroot (only used when dir is provided). Only match stat files named <outroot>.stat
##' @param stats A list of statistics to plot (i.e., "likelihood", "recombs", "prior", etc).
##' If NULL, then plot all statistics found in the stat files
##' @param colors A list of colors, one per file, to use in the plots
##' @param xlim The range of iterations to plot. If NULL, use entire range
##' @param callPar If TRUE, call the par() function to create a new graphics object with
##' a plot for each statistic on the same page
##' @param ncol (Only used if callPar=TRUE), arrange the plots in this many columns
##' @param addLegend If TRUE, add an additional plot that contains a legend mapping files to colors
##' @param ... Arguments passed to plotStats()
##' @return invisibly returns a data.frame detailing the files and colors used
##' @export
plotTraces <- function(files=NULL, dir=NULL, outroot=NULL, stats=NULL,
                       colors=if(length(files)==1) { "black"} else { rainbow(length(files))},
                       xlim=NULL, callPar=TRUE, ncol=2, addLegend=(length(files) != 1),
                       ...) {
    if (is.null(files)) {
        if (is.null(dir)) stop("either files or dir must be provided")
        if (!is.null(outroot)) pattern <- sprintf("%s.stats", outroot) else pattern <- ".*.stats"
        files <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
        if (length(files)==0) stop("No stat files found in ", dir)
        if (missing(colors)) colors <- rainbow(length(files))
    }
    if (length(files) != length(colors)) stop("number of files/colors do not match")
    statsGiven <- !is.null(stats)
    xlimGiven <- !is.null(xlim)
    x <- list()
    for (i in 1:length(files)) {
        currx <- readStatsFile(files[i])
        if (!xlimGiven) {
            xlim <- range(c(xlim, currx$iter))
        } else {
            currx <- currx[currx$iter >= xlim[1] & currx$iter <= xlim[2],]
        }
        if (!statsGiven)
            stats <- unique(c(stats, names(currx[,3:ncol(currx)])))
        x[[i]] <- currx
    }
    if (callPar) {
        nrow <- ceiling((length(stats)+addLegend)/ncol)
        cat(nrow, ncol, "\n")
        par(mfrow=c(nrow, ncol), mar=c(4,4,1,1), mgp=c(2,0.75,0))
    }
    cat(xlim, "\n")
    cat(class(xlim), "\n")
    cat(colors, "\n")
    for (stat in stats) {
        cat("plotting stat", stat, "\n")
        ylim <- NULL
        for (i in 1:length(files))
            if (is.element(stat, names(x[[i]])))
                ylim <- range(c(ylim, x[[i]][,stat]))
        add <- FALSE
        for (i in 1:length(files)) {
            if (is.element(stat, names(x[[i]]))) {
                plotStats(x[[i]], stat, add=add, col=colors[i], xlim=xlim, ylim=ylim, ...)
                add <- TRUE
            }
        }
    }
    legendData <- data.frame(file=files, color=colors)
    if (addLegend) {
        tmplist <- strsplit(files, "/")
        maindir <- NULL
        while(1) {
            tmp <- unique(sapply(tmplist, function(x) {x[1]}))
            if (length(tmp) == 1) {
                maindir <- c(maindir, tmp)
                tmplist <- lapply(tmplist, function(x) {x[-1]})
            } else break
        }
        shortFiles <- sapply(tmplist, function(x) {paste(x, sep="", collapse="/")})
        plot(c(0),c(0),type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
        legend(x="topleft", shortFiles, pch=1, col=colors)
        if (!is.null(maindir))
            mtext(sprintf("base dir = %s", paste(maindir, sep="", collapse="/")))
    }
    invisible(legendData)
}

