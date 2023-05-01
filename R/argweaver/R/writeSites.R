
##' Write an ARGweaver sites file
##' @param sites Sites object
##' @param file Sites file (may be gzipped. may NOT be streamed)
##' @return a list with two elements: one describing the region, and the
##' other a data frame with a column for each haploid genome and a row
##' for each variant site. A third element "probs" will be added if readProbs=TRUE
##' @export
writeSites <- function(sites, filename) {
    if (grepl(".gz$", filename)) {
        outfile <- gzfile(filename, "w")
    } else outfile <- file(filename, "w")
    cat(c("#NAMES", names(sites$sites)), file=outfile, sep="\t")
    cat("\n", file=outfile)
    cat(c("#REGION", as.character(sites$region$chrom), format(c(sites$region$start, sites$region$end), scientific=FALSE)),
        file=outfile, sep="\t")
    cat("\n", file=outfile)
    sitePatterns <- apply(sites$sites, 1, paste, collapse="")
    df <- data.frame(pos=sites$pos, pattern=sitePatterns)
    if (!is.null(sites$probs)) {
        df <- cbind(df, sites$probs)
    }
    write.table(format.data.frame(df, scientific=FALSE, digits=3, trim=TRUE),
                file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    close(outfile)
    invisible(NULL)
}



