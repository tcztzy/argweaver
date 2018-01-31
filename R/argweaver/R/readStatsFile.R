
##' Read arg-sample stats file
##'
##' This function will read a <outroot>.stats file written by arg-sample.
##' @param file A character string, either the name of the stats file or
##' the directory where the stat file can be found.
##' If a directory, there must be only one
##' stat file inside or else an error will occur
##' @param burnin If not NULL, a single numeric value. Any statistics from
##' iterations < burnin will be discarded
##' @param thinVal Remove all iters for which iteration mod thinVal != 0
##' @return A data frame containing the data from the stats file. Only
##' MCMC iterations of type "resample" will be retained. In addition to
##' the original columns from the stats file, several statistics may be
##' added to the data frame as follows:
##' "joint2" is the sum of prior2 and likelihood and is an alternative
##' posterior likelihood statistic that is preferred over "joint",
##' which is the sum of prior and likelihood. Whereas "prior" is the
##' prior probability of the sampled ARG under the model, "prior2"
##' integrates over all possible recombination events that could produce
##' the same series of local trees.
##' "total_recomb" is the sum of "recombs" and "invis_recombs"
##' @export
readStatsFile <- function(file, burnin=NULL, thinVal=1) {
    if (length(file) != 1)
        stop("readStatsFile can currently only read single file")
    if (!file.exists(file))
        stop(file, " does not exist")
    if (!is.null(burnin)) {
        if (length(burnin) != 1)
            stop("burnin should be length 1")
        if (class(burnin) != "numeric" && class(burnin) != "integer")
            stop("burnin should be a numeric value")
    }
    if (file.info(file)$isdir) {
        possibleFiles <- list.files(file, pattern=".*.stats",
                                    recursive=TRUE, full.names=TRUE)
        if (length(possibleFiles) == 0)
            stop(file, "is a directory but no stats file found within")
        if (length(possibleFiles) != 1)
            stop(file, " is a directory and multiple stats files found within")
        file <- possibleFiles
    }
    x <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    x <- x[x$stage=="resample",]
    if (!is.null(burnin))
        x <- x[x$iter >= burnin,]
    x <- x[x$iter%%thinVal == 0,]
    if (is.element("prior2", names(x)))
        x$joint2 <- x$prior2 + x$likelihood
    if (is.element("recombs", names(x)) &&
        is.element("invis_recombs", names(x)))
        x$total_recombs <- x$recombs + x$invis_recombs
    if (substr(file, nchar(file)-5, nchar(file))==".stats") {
        fixFile <- sprintf("%s.updated_likelihoods.bed", substr(file, 1, nchar(file)-6))
        if (!file.exists(fixFile)) return(x)
        newx <- read.table(fixFile, header=TRUE, stringsAsFactors=FALSE)
        if (is.element("rep", names(newx))) names(newx)[which(names(newx)=="rep")] <- "iter"
        x <- x[is.element(x$iter, newx$iter),]
        if (nrow(x) != nrow(newx) || sum(x$iter != newx$iter) != 0) stop("Error fixing stats")
        x$likelihood <- newx$likelihood
        x$newPrior <- newx$prior
        x$newPrior2 <- newx$prior2
        x$newRecombs <- newx$recombs
        x$newNoncompats <- newx$noncompats
        x$joint <- x$prior + x$likelihood
        x$joint2 <- x$prior2 + x$likelihood
    }
    x
}
