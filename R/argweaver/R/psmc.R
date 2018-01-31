##' Make a popsize file suitable for arg-sample from PSMC output
##' @param psmcFile File with PSMC output
##' @param ntimes Number of time points in arg-sample model
##' @param delta arg-sample delta parameter
##' @param maxTime maximum time point used in arg-sample
##' @param mu mutation rate in mutations per base pair per generation
##' @param s The bin size used in the fq2psmcfa command; the default in that program is 100
##' @param sampleAge The age of the sample, in generations (For ancient samples)
##' @param outfile If not NULL, path to write output popsize file. This file is suitable
##' for use with the "--popsize-file" option in arg-sample.
##' @return a data frame with time and popsize (invisibly if outfile is provided)
##' @export
##' @import utils
psmcToPopsize <- function(psmcFile, ntimes=20, delta=0.01, maxTime=200000,
                          mu=1.4e-8, s=100,
                          sampleAge=0, outfile=NULL) {
    x <- read.table(psmcFile, header=FALSE, fill=TRUE, stringsAsFactors=FALSE,
                    col.names=sprintf("col%i", 1:100))
    theta <- as.numeric(tail(x[x$col1=="TR",2], 1))
    rs <- x[x$col1=="RS",2:7]
    lastRs <- tail(which(rs$col2 == "0"), 1)
    rs <- rs[lastRs:nrow(rs),]
    names(rs) <- c("k", "t", "lambda", "pi", "ahet", "ahom")
    for (i in 1:ncol(rs)) rs[,i] <- as.numeric(rs[,i])
    
    n0 <- theta / (4*mu) / s

    rs$minTime <- 2*n0*rs$t
    rs$minTime[-1] <- rs$minTime[-1] + sampleAge
    rs$maxTime <- c(rs$minTime[-1], Inf)
    rs$popsize <- n0 * rs$lambda

    ## Now we have psmc's estimates. Need to convert them to ARGweaver
    times <- getTimes(n=ntimes, maxTime=maxTime, delta=delta)
    halfTimes <- getTimes(n=2*ntimes-1, maxTime=maxTime, delta=delta)

    results <- NULL
    lastTime <- 0
    for (i in unique(c(seq(from=2, by=2, to=length(halfTimes)), length(halfTimes)))) {
        minTime <- lastTime
        maxTime <- halfTimes[i]
        lastTime <- halfTimes[i]
        tmp <- rs[rs$minTime < maxTime & rs$maxTime > minTime,]
        tmp[1,"minTime"] <- minTime
        tmp[nrow(tmp),"maxTime"] <- maxTime
        tmp$tlen <- tmp$maxTime - tmp$minTime

        popsize <- sum(tmp$tlen)/sum(tmp$tlen / tmp$popsize)
        tmp2 <- data.frame(t=maxTime, popsize=popsize)
        if (is.null(results)) results <- tmp2 else results <- rbind(results, tmp2)
    }
    if (is.null(outfile)) return(results)
    write.table(results, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    invisible(results)
}


##' For multiple population model in arg-sample, create popsize file using
##' one representative PSMC result file per population
##' @param psmcFiles A list with length = number of populations. Each element should
##' be a character string giving the file path of PSMC result for a population,
##' or can be a single numeric value for constant size population (giving diploid
##' population size)
##' @param ages (Used for ancient samples) A numeric vector of same length as psmcFiles.
##' Each value represents the age in generations of the sample used in PSMC.
##' @inheritParams psmcToPopsize
##' @return A data frame with columns: population, time, popsizse (invisibly if outfile is
##' not NUL)
##' @note Note that this function is for multiple-population model of arg-sample, which is still
##' experimental. Most users will want a single vector of population sizes, using the
##' function psmcToPopsize function instead.
##' @export
multiPopsizeFileFromPsmc <- function(psmcFiles, ages=NULL, ntimes=20,
                                     delta=0.01, maxTime=1e6,
                                     mu=1.4e-8, s=100, outfile=NULL) {
    rv <- NULL
    if (!is.list(psmcFiles))
        stop("first argument should be a list of length = number of populations")
    for (i in 1:length(psmcFiles)) {
        if (is.numeric(psmcFiles[[i]])) {
            results <- data.frame(pop=i-1, t=maxTime, popsize=psmcFiles[[i]])
        } else {
            if (!is.character(psmcFiles[[i]]))
                stop("each element of psmcFiles should be a popsize or a path to psmc resutls")
            if (is.null(ages)) sampleAge <- 0 else sampleAge <- ages[i]
            results <- psmcToPopsize(psmcFiles[[i]],
                                 ntimes=ntimes, delta=delta, maxTime=maxTime,
                                 mu=mu, s=s, sampleAge=sampleAge)
            results <- data.frame(pop=i-1, t=results$t, popsize=results$popsize)
        }
        if (is.null(rv)) rv <- results else rv <- rbind(rv, results)
    }
    if (is.null(outfile)) return(results)
    write.table(results, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    invisible(results)
}


