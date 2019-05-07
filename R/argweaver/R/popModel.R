
meanPopPerBranch <- function(idx, edge, mod) {
    t1 <- getDiscreteTime(edge[idx,"depth"], mod)
    t2 <- getDiscreteTime(edge[idx,"depth"]+edge[idx,"length"], mod)
    if (t1 == t2) return(c(0,0))
    val <- 0
    path <- edge$path[idx]
    for (time in t1:(t2-1)) {
        pop <- mod$paths[path+1,time]
        val <- val + pop * (mod$times[time+1] - mod$times[time])
    }
    c(val, edge[idx,"length"])
}

meanPopRec <- function(idx, edge, mod) {
    node <- edge[idx,"dec"]
    anc <- which(edge$anc==node)
    thisVal <- meanPopPerBranch(idx, edge, mod)
    if (length(anc)==0) {  #leaf
        return(thisVal)
    }
    if (length(anc) != 2) stop("Error")
    return(thisVal + meanPopRec(anc[1], edge, mod) + meanPopRec(anc[2], edge, mod))
}

meanPop <- function(idx, edge, mod) {
    val <- meanPopRec(idx, edge, mod)
    if (val[2] == 0) return(0)
    val[1]/val[2]
}



##' scalePopModel: scale times in a population model
##' @param mod A popmodel (usually read by readPopModel)
##' @param scale Scaling factor
##' @return a new popmodel with scaled times
##' @export
scalePopModel <- function(mod, scale) {
    mod$times <- mod$times*scale
    mod$div$time <- mod$div$time*scale
    if (!is.null(mod$mig))
        mod$mig$time <- mod$mig$time*scale
    mod
}


drawPopModel <- function(mod, add=FALSE, col=rgb(0,0.8,0.3,0.15), logscale=FALSE, popwidth=NULL, 
                         betweenPopWidth = 0.2, popnames=NULL, divwidth=0.01, migCol=rgb(1,0,0,0.5),
                         migFactor=0.5,
                         minMig=0.01, drawTimes=TRUE, ylim=NULL,
                         timescale=1, ylab="Time", xlab="Population") {
    if (timescale != 1) mod <- scalePopModel(mod, timescale)
    if (is.null(popwidth)) popwidth <- rep(0.8, mod$npop)
    npop <- mod$npop
    if (is.null(popnames)) popnames <- 1:npop
    if (!add) {
        usetimes <- mod$times
        if (logscale) usetimes[usetimes==0] <- 1
        if (is.null(ylim)) ylim <- range(usetimes)
        xlim <- c(0, sum(popwidth) + betweenPopWidth*(mod$npop+1))
        plot(xlim,ylim, type="n",
             xlim=xlim, ylim=ylim,
             xlab=xlab, ylab=ylab,
             log=ifelse(logscale, "y", ""), xaxt="n", yaxs="i")
    }

    popCoords <- NULL
    xend <- 0
    for (pop in 1:mod$npop) {
        endTime <- mod$div[mod$div[,"pop1"]==pop-1,"time"]
        if (length(endTime)==0) {
            endTime <- max(mod$times)
        }
        xstart <- xend + betweenPopWidth
        xend <- xstart + popwidth[pop]
        coords <- data.frame(pop=pop-1, x0=xstart, x1=xend, y0=ifelse(logscale, 1, 0), y1=endTime)
        popCoords <- rbind(popCoords, coords)
        rect(xleft=coords$x0, xright=coords$x1, ybottom=coords$y0, ytop=coords$y1, col=col, border=NA)
##        mtext(text=popnames[pop], side=1, at=(xstart+xend)/2, line=0.5)
    }
    
    ## draw div lines
    mod$div$y0=-1
    mod$div$y1=-1
    divwidth <- divwidth * (par("usr")[4]-par("usr")[3])    
    for (i in 1:nrow(mod$div)) {
        pop1 <- mod$div[i,"pop1"]
        pop2 <- mod$div[i,"pop2"]
        ## ensure pop1 < pop2
        if (pop1 > pop2) {
            tmp <- pop2
            pop2 <- pop1
            pop1 <- tmp
        }
        if (logscale) {
            mod$div$y1[i] <- 10^(log10(mod$div[i,"time"]))
            mod$div$y0[i] <- 10^(log10(mod$div[i,"time"])-divwidth)
        } else {
            mod$div$y1[i] <- mod$div[i,"time"]
            mod$div$y0[i] <- mod$div[i,"time"]-divwidth
        }
        rect(xleft=popCoords[popCoords$pop==pop1,"x1"],
             xright=popCoords[popCoords$pop==pop2,"x0"],
             ytop=mod$div$y1[i], ybottom=mod$div$y0[i], col=col, border=NA)
    }
    if (!is.null(mod$mig)) {
        mod$mig$y0 <- -1
        mod$mig$y1 <- -1
        for (i in 1:nrow(mod$mig)) {
            pop1 <- mod$mig[i,"pop1"]
            pop2 <- mod$mig[i,"pop2"]
            if (pop1 < pop2) {
                xFrom <- popCoords[popCoords$pop==pop1,"x1"]
                xTo <- popCoords[popCoords$pop==pop2,"x0"]
            } else {
                xFrom <- popCoords[popCoords$pop==pop1,"x0"]
                xTo <- popCoords[popCoords$pop==pop2,"x1"]
            }
            if (logscale) {
                mod$mig$y1[i] <- 10^(log10(mod$mig[i,"time"])+divwidth/2)
                mod$mig$y0[i] <- 10^(log10(mod$mig[i,"time"])-divwidth/2)
            } else {
                mod$mig$y1[i] <- mod$mig[i,"time"]+divwidth/2
                mod$mig$y0[i] <- mod$mig[i,"time"]-divwidth/2
            }
            arrowLen <- max(mod$mig[i,"prob"]*migFactor, minMig)
            arrows(x0=xFrom,x1=xTo,y0=mod$mig[i,"time"], col=migCol,
                   length=par("pin")[1]*arrowLen)
        }
    }
    if (drawTimes) abline(h=mod$times[-1], lty=3, col="gray")
    ## TODO: add mig bands and connect pops at div times
    mod$popCoords <- popCoords
    invisible(scalePopModel(mod, 1/timescale))
}


##' Get population paths from ARGweaver log file
##' @param logfile The log file name
##' @param ntimes The number of discrete times
##' @return A data frame with one row for each path; each column is a time interval
##' @export
getPathsFromLogFile <- function(logfile, ntimes) {
    x <- scan(logfile, what="character", nmax=10000, quiet=TRUE)
    w <- which(x=="numpath")
    if (length(w) != 1) stop("error1")
    numpath <- as.numeric(x[w+2])
    if (numpath < 1) stop("error2")
    paths <- data.frame(matrix(nrow = numpath, ncol = ntimes))
    rownames(paths) <- sprintf("path%i", (1:numpath)-1)
    colnames(paths) <- sprintf("time%i", 1:ntimes)
    for (i in 1:numpath) {
        pathname <- sprintf("path%i", i-1)
        w <- which(x == pathname)
        if (length(w) != 1) stop("error3")
        path <- x[(w+2):(w+2+ntimes-1)]
        path[1] <- substr(path[1], 2, nchar(path[1]))
        for (j in 1:ntimes) path[j] <- substr(path[j], 1, nchar(path[j])-1)
        path <- as.numeric(path)
        paths[pathname,] <- path
    }
    paths
}





##' Read in a population tree
##' @param popfile the file passed to arg-sample wtih --pop-tree-file
##' @return A list with npop, divergence events (divs) and migration events (migs)
##' @export
readPopFile <- function(popfile) {
    x <- read.table(popfile, fill=TRUE, stringsAsFactors=FALSE, col.names=c("event", "time", "pop1", "pop2", "prob"))
    if (x[1,1] != "npop") {
        stop("Error1 reading popfile")
    }
    npop <- as.numeric(x[1,2])
    divs <- x[x[,1]=="div",2:4]
    if (nrow(divs) != npop-1) stop("Error2 reading popfile")
    names(divs) <- c("time", "pop1", "pop2")
    divs$time <- as.numeric(divs$time)
    if (ncol(x) >= 5) {
        migs <- x[x[,1]=="mig",2:5]
        names(migs) <- c("time", "pop1", "pop2", "prob")
        migs$time <- as.numeric(migs$time)
        migs$prob <- as.numeric(migs$prob)
        if (nrow(migs)==0) migs <- NULL
    } else migs <- NULL
    list(npop=npop, div=divs, mig=migs)
}


##' Read in a population model file
##' @param logfile The log file output by arg-sample
##' @param popfile The population file given to arg-sample with the argument
##' --pop-tree-file
##' @return A population object suitable to pass to plotTree function
##' @export
readPopModel <- function(logfile, popfile) {
    mod <- readPopFile(popfile)
    x <- scan(logfile, what="character", nmax=10000, quiet=TRUE)
    w <- which(x=="Using")
    for (i in w) {
        if (x[i+1] != "time") next
        tNew <- as.numeric(x[i+2])
        tOld <- as.numeric(x[i+5])
        event <- x[i+7]
        if (event == "div") {
            f <- mod$div$time == tOld
            if (sum(f) > 0) {
                ## if there are multiple events at the same time they may already have been fixed
                mod$div[mod$div$time==tOld,"time"] <- tNew
            }
        }
        if (event == "mig") {
            f <- mod$mig$time == tOld
            ## if there are multiple events at the same time they may already have been fixed
            if (sum(f) > 0) {
                mod$mig[mod$mig$time==tOld,"time"] <- tNew
            }
        }
    }
    mod$times <- getTimesFromLogFile(logfile)
    mod$ntimes <- length(mod$times)
    mod$paths <- getPathsFromLogFile(logfile, mod$ntimes)
    mod$npaths <- nrow(mod$paths)
    mod
}
