
##' Read arg-sample summary file
##'
##' This function will read a file written by arg-summarize
##' @param file A character string, either the name of the stats file or
##' the directory where the stat file can be found.
##' If a directory, there must be only one
##' stat file inside or else an error will occur
##' @param burnin If not NULL, a single numeric value. Any statistics from
##' iterations < burnin will be discarded. Will be ignored if there is not
##' a column named "MCMC_sample"
##' @return A data frame containing the arg-summarize results
##' @export
readArgSummary <- function(file, burnin=NULL) {
    if (!is.element("pipe", class(file))) {
        if (length(file) != 1)
            stop("readSummary can currently only read single file")
        if (!file.exists(file))
            stop(file, " does not exist")
    }
    x <- read.table(file, header=TRUE, skip=2, comment.char="")
    names(x)[1] <- "chrom"
    if (!is.null(burnin)) {
        if (length(burnin) != 1)
            stop("burnin should be length 1")
        if (class(burnin) != "numeric" && class(burnin) != "integer")
            stop("burnin should be a numeric value")
        if (!is.element("MCMC_rep", names(x)))
            warning("Ignoring burnin argument; MCMC rep not in file")
        x <- x[x$MCMC_rep >= burnin,]
    }
    x
}


##' Choose columns from summary file
##'
##' Given an ARGweaver summary file,
##' return columns of interest
##' @param x A data frame read by readSummary()
##' @param inds A vector of character strings containing individual
##' (or haplotype names) relevant to desired columns
##' @param prefix Choose only columns with this prefix
##' @param suffix Choose only columns with this suffix
##' @param or If TRUE, return columns that mention any of the individuals.
##' If FALSE, return columns that mention all the individuals
##' @return A list of column names from x that meet the criteria
##' @note This may not work as desired if some individual names are subsets of others
##' @export
getCols <- function(x, inds, prefix="", suffix="", or=FALSE) {
    x <- names(x)
    if (nchar(prefix) > 0)
        x <- grep(sprintf("^%s", prefix), x, value=TRUE)
    if (nchar(suffix) > 0)
        x <- grep(sprintf("%s$", suffix), x, value=TRUE)
    if (length(x) == 0) return(x)
    if (!or) {
        for (ind in inds)
            x <- grep(ind, x, value=TRUE)
    } else {
        keep <- rep(FALSE, length(x))
        for (ind in inds)
            keep <- keep | grepl(ind, x)
        x <- inds[keep]
    }
    x
}

