
##' Plot a sites file
##' @param x A sites file object as read by readSites().
##' @param file (Alternative to x) Sites file (may be gzipped. may NOT be streamed)
##' @param start The start coordinate of the plot. NULL implies start from beginning of region
##' @param end The end coordinate of the plot. NULL implies go to end of region.
##' @param anc Ancestral sample name. If NULL, allele1 is major allele and allele2 is minor allele
##' @param singletons If TRUE, include singletons in plot
##' @param allele1Color The color(s) for allele1 (will be recycled across sites)
##' @param allele2Color The color(s) for allele2 (will be recycled across sites)
##' @param missingColor The color for N allleles
### @param alleleColors A list giving color for each allele. alleleColors$anc gives
### color for allele matching ancestral sample
##' @param stretch A horizontal expansion factor for sites (so that sparse sites are
##' visible in a large region)
##' @param indOrder A character vector giving order to list haplotypes (if NULL use
##' order used in sites object)
##' @param useCoords if FALSE, just plot all variant sites with equal spacing
##' @param xlab label for x-axis. Default is "Coordinate" if useCoords=TRUE
##' or "Sites index" if useCoords=FALSE
##' @param textColor a list saying which color to use for which haplotype label
##' @param textSize the magnification size of haplotype labelscat 
##' @param ... Other arguments passed to plot function
##' @export
plotSites <- function(x=NULL, file=NULL, start=NULL, end=NULL,
#                      baseColor=list(N=rgb(0,0,0,0.2), A="red",
                                        #                                     C="green", G="blue", T="purple"),
                      anc=NULL,
                      singletons=FALSE,
#                      alleleColors=list(A="red", C="green", G="blue", T="purple", N=rgb(0,0,0,0.1), anc="black"),
                      allele1Color="black",
                      allele2Color=c("red", "orange", "yellow", "green", "turquoise", "purple"),
                      missingColor=rgb(0,0,0,0.1),
                      stretch=1, indOrder=NULL, textColor=NULL,
                      useCoords=TRUE, xlab=NULL, textSize=1, ...) {
    if (is.null(x)) {
        if (is.null(file)) stop("Must provide either x or file\n")
        x <- readSites(file)
    }
    if (is.null(start)) start <- x$region$start
    if (is.null(end)) end <- x$region$end

    f <- x$pos >= start & x$pos <= end
    if (sum(f) == 0) stop("No sites in start=",start," end=",end)
    x$pos <- x$pos[f]
    x$sites <- x$sites[f,]

    numind <- ncol(x$sites)
    if (!is.null(indOrder)) {
        ypos <- numeric(numind)
        for (i in 1:numind) {
            tmp <- which(indOrder == names(x$sites)[i])
            if (length(tmp) == 0) ypos[i] <- NA else ypos[i] <- tmp
        }
        x$sites <- x$sites[,!is.na(ypos)]
        ypos <- ypos[!is.na(ypos)]
        numind <- length(ypos)
    } else {
        ypos <- 1:numind
    }
    count <- list()
    alleles <- c("A", "C", "G", "T")
    if (is.null(anc)) {
        tmpsites <- x$sites
    } else {
        tmpsites <- x$sites[,names(x$sites) != anc]
    }
    for (allele in alleles)
        count[[allele]] <- apply(tmpsites, 1, function(x) {sum(as.character(x)==allele)})
    if (singletons) {minCount <- 1} else {minCount <- 2}
    f <- (( (count[["A"]] >= minCount ) +
            (count[["C"]] >= minCount ) +
            (count[["G"]] >= minCount ) +
            (count[["T"]] >= minCount )) > 1)
    x$pos <- x$pos[f]
    x$sites <- x$sites[f,]
    if (useCoords) {
        xmin <- x$pos-stretch/2
        xmax <- x$pos+stretch/2
    } else {
        xmin <- seq(from=0.5, to=nrow(x$sites)-0.5, by=1)
        xmax <- xmin+1
    }
    for (allele in alleles) count[[allele]] <- count[[allele]][f]
    majAllele <- mapply(function(a,c,g,t) {x <- c(a,c,g,t)
        names(x) <- alleles
        names(x)[which.max(x)]}, count[["A"]], count[["C"]], count[["G"]], count[["T"]])
    if (!is.null(anc)) {
        w <- which(names(x$sites) == anc)
        if (length(w) != 1) stop("Error finding ancestral sequence ", anc)
        majAllele <- ifelse(x$sites[,w] == "N", majAllele, as.character(x$sites[,w]))
    }
    par(mar=c(4,6,3,1), mgp=c(2,1,0), xaxs="i", yaxs="i")

    ylim <- c(min(ypos,na.rm=TRUE)-0.5, max(ypos,na.rm=TRUE)+0.5)
    if (is.null(xlab)) {
        if(useCoords) {xlab <- "Coordinate"
        } else {
            xlab <- "Site index"
        }
    }
    plot(c(0),c(0), type="n", yaxt="n", xlim=c(min(xmin), max(xmax)),
         ylim=ylim, xlab=xlab, ylab="", ...)
    allele2Color <- rep(allele2Color, length.out=nrow(x$sites))
    allele1Color <- rep(allele1Color, length.out=nrow(x$sites))
    for (i in 1:numind) {
        if (is.na(ypos[i])) next
        rect(xleft=xmin, xright=xmax,
             ybottom=ypos[i]-0.5, ytop=ypos[i]+0.5, border=NA,
             col=ifelse(x$sites[,i]==majAllele,allele1Color,ifelse(x$sites[,i]=="N", missingColor, allele2Color)))
        if (is.null(textColor) || is.null(textColor[[names(x$sites)[i]]])) {
            textcol <- "black"
        } else textcol <- textColor[[names(x$sites)[i]]]
        mtext(names(x$sites)[i], side=2, at=ypos[i], las=2, col=textcol, cex=textSize)
    }
    abline(h=seq(from=min(ypos, na.rm=TRUE)-0.5, to=max(ypos,na.rm=TRUE)+0.5, by=2), lty=3, col="gray")
    invisible(x)
}


##' call plotTreeFromBed at the sites in a site file; color individual nodes by allele
##' @param bedFile bed.gz file from smc2bed with trees
##' @param sites Sites object, as returned by readSites(). Will plot one tree for each site.
##' @param pos If not NULL, a vector of integers giving which sites to plot.
##' If NULL, plot all sites
##' @param start start position
##' @param end end position
##' @param doSingletons whether to plot singletons
##' @param colors Color to use for each allele in leaf of tree
##' @param ... Additional arguments passed to plotTreesFromBed
##' @export
plotTreesAtSites <- function(bedFile, sites, pos=NULL, start=NULL, end=NULL, doSingletons=FALSE,
                             colors=list(A="red", C="green", G="turquoise", T="purple",N="gray"),
                             ...) {
    if (is.null(pos)) pos <- sites$pos
    if (!is.null(start)) pos <- pos[pos >=start]
    if (!is.null(end)) pos <- pos[pos <= end]
    numskip <- 0
    for (i in 1:length(pos)) {
        f <- which(sites$pos == pos[i])
        if (length(f) != 1) {
            numskip <- numskip+1
            next
        }
        tmp <- table(as.character(sites$sites[f,]))
        tmp <- tmp[names(tmp)!="N"]
        if (length(tmp) == 1) next
        if (!doSingletons) {
            if (sum(tmp > 1) <= 1) next
        }
        currCol <- list()
        currCol[names(sites$sites)] <- colors[as.character(unlist(sites$sites[f,]))]
        cat("plotting pos ", pos[i], "\n")
        plotTreesFromBed(bedFile, col=currCol, start=pos[i], end=pos[i], treeInfo=pos[i], ...)
    }
    if (numskip > 0) warning("Skipped ", numskip, " sites that were not found")
    invisible(NULL)
}
