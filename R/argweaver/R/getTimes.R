
##' Compute log-distributed discrete times used by arg-sample
##'
##' This function returns a vector of times which are spaced
##' on a log scale according to the parameter delta.
##' @param n The number of time points between 0 and maxTime, inclusive
##' @param maxTime The maximum time
##' @param delta Spacing parameter. Times get denser towards t=0
##' as delta gets small
##' @return A vector of times of length n
##' @export
getTimes <- function(n=20, maxTime=200000, delta=0.01) {
  (exp((0:(n-1))/(n-1) * log(1 + delta*maxTime)) -1)/delta
}




##' Determine argweaver log file starting with file name
##' @param filename Beginning of log filename
##' @return name of ARGweaver log file if there is one good candidate, otherwise error
##' @export
getLogFile <- function(filename) {
    dir <- dirname(filename)
    candidates <- list.files(path=dir, full.names=TRUE, pattern=".*.log")
    if (length(candidates)==0) return(NULL)
    if (length(candidates)==1) return(candidates)
    bn <- basename(filename)
    dotPos <- gregexpr('.', bn, fixed=TRUE)[[1]]
    if (length(dotPos)==1 && dotPos == -1) {
        stop("cannot determine log filename")
    }
    candidates <- NULL
    for (i in 1:length(dotPos)) {
        tmp <- sprintf("%s/%s.log", dir, substr(bn, 1, dotPos[i]-1))
        if (file.exists(tmp))
            candidates <- c(candidates, tmp)
    }
    if (length(candidates) > 1)
        stop("got multiple log file candidates\n")
    if (length(candidates) == 0)
        stop("cannot determine log file")
    candidates
}


##' Get times from an argweaver log file
##' @param logfile the .log file output by ARGweaver
##' @return Numeric list of discrete times
##' @export
getTimesFromLogFile <- function(logfile) {
    if (substr(logfile, nchar(logfile)-3, nchar(logfile)) != ".log")
        logfile <- getLogFile(logfile)
    x <- scan(logfile, what="character", nmax=10000, quiet=TRUE)
    w <- which(x=="times")[1]
    if (x[w+1] != "=") stop("error reading logfile")
    w <- w+2
    str <- substr(x[w], start=2, stop=nchar(x[w])-1)
    as.numeric(strsplit(str, ",")[[1]])
}


