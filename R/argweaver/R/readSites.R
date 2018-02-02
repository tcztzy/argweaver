
##' Read an ARGweaver sites file
##' @param file Sites file (may be gzipped. may NOT be streamed)
##' @return a list with two elements: one describing the region, and the
##' other a data frame with a column for each haploid genome and a row
##' for each variant site
##' @export
readSites <- function(file) {
    names <- scan(file, nlines=1, what=character(), quiet=TRUE)[-1]
    region <- scan(file, nlines=2, what=character(), quiet=TRUE)[-(1:(length(names)+1))]
    if (region[1] != "REGION" && region[1] != "#REGION")
        stop("Expected string REGION at beginning of SITES file\n")
    if (length(region) != 4)
        stop("Format error in SITES file\n")
    rv <- list()
    rv$region <- data.frame(chrom=region[2], start=as.numeric(region[3]), end=as.numeric(region[4]))
    x <- read.table(file, skip=2, header=FALSE, stringsAsFactors=FALSE)[,1:2]
    rv$pos <- x[,1]
    tmp <- unlist(strsplit(x[,2], ""))
    for (i in 1:length(names))
        x <- cbind(x, tmp[seq(from=i, to=length(tmp), by=length(names))])
    x <- x[,-c(1:2)]
    names(x) <- c(names)
    rv$sites <- x
    rv
}
