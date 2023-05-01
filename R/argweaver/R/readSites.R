
##' Read an ARGweaver sites file
##' @param file Sites file (may be gzipped. may NOT be streamed)
##' @param readProbs If TRUE, return genotype probabilities as well (if they are in encoded in the file)
##' @return a list with two elements: one describing the region, and the
##' other a data frame with a column for each haploid genome and a row
##' for each variant site. A third element "probs" will be added if readProbs=TRUE
##' @export
readSites <- function(file, readProbs=FALSE) {
    names <- scan(file, nlines=1, what=character(), quiet=TRUE)[-1]
    region <- scan(file, nlines=2, what=character(), quiet=TRUE)[-(1:(length(names)+1))]
    if (region[1] != "REGION" && region[1] != "#REGION")
        stop("Expected string REGION at beginning of SITES file\n")
    if (length(region) != 4)
        stop("Format error in SITES file\n")
    rv <- list()
    rv$region <- data.frame(chrom=region[2], start=as.numeric(region[3]), end=as.numeric(region[4]))
    cols <- 1:2
    if (readProbs) cols <- 1:(4*length(names)+2)
    x <- read.table(file, skip=2, header=FALSE, stringsAsFactors=FALSE)[,cols]
    if (ncol(x) == 2) readProbs <- FALSE
    if (readProbs) {
        probs <- list()
        idx <- 3
        for (ch in c("A", "C", "G", "T")) {
            tmp <- x[,seq(from=idx, to=ncol(x), by=4)]
            names(tmp) <- names
            probs[[sprintf("prob%s", ch)]] <- tmp
            idx <- idx+1
        }
        x <- x[,1:2]
    }
    rv$pos <- x[,1]
    tmp <- unlist(strsplit(x[,2], ""))
    for (i in 1:length(names))
        x <- cbind(x, tmp[seq(from=i, to=length(tmp), by=length(names))], stringsAsFactors=FALSE)
    x <- x[,-c(1:2)]
    names(x) <- c(names)
    rv$sites <- x
    if (readProbs) {
        rv$probs <- probs
    }
    rv
}



#numAllele <- function(x) {
#    
#}


##' Subset a sites object
##' @param x a sites object (read in by readSites)
##' @param keep A list of haplotypes to keep
##' @param prune A list of haplotypes to prune
##' @param region A numeric vector of length 2 giving minimum and maximum coordinates
##' (1-based)
##' @return A sites object subsetted by haplotypes and/or region
##' @export
subSites <- function(x, keep=NULL, prune=NULL, region=NULL) {
    if (!is.null(keep)) {
        f <- is.element(names(x$sites), keep)
    } else f <- rep(TRUE, length(names(x$sites)))
    if (!is.null(prune))
        f <- f & !is.element(names(x$sites), prune)
    x$sites <- x$sites[,f]
    numUnique <- apply(x$sites, 1, function(x) {length(unique(x))})
    f <- numUnique > 1 | x$sites[,1]=="N"
    x$sites <- x$sites[f,]
    x$pos <- x$pos[f]
    if (!is.null(region)) {
        if (length(region) != 2) stop("region should be numeric vector of length 2")
        f <- x$pos >= region[1] & x$pos <= region[2]
        x$pos <- x$pos[f]
        x$sites <- x$sites[f,]
        x$region$start <- region[1]
        x$region$end <- region[2]
    }
    x
}
