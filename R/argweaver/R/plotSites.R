
##' Plot a sites file
##' @param x A sites file object as read by readSites().
##' @param file (Alternative to x) Sites file (may be gzipped. may NOT be streamed)
##' @param start The start coordinate of the plot. NULL implies start from beginning of region
##' @param end The end coordinate of the plot. NULL implies go to end of region.
##' @param anc Ancestral sample name. If NULL, allele1 is major allele and allele2 is minor allele
##' @param allele1Color The color(s) for allele1 (will be recycled across sites)
##' @param allele2Color The color(s) for allele2 (will be recycled across sites)
##' @param missingColor The color for N allleles
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
    for (allele in alleles)
        count[[allele]] <- apply(x$sites, 1, function(x) {sum(as.character(x)==allele)})
    f <- (( (count[["A"]] > 1 ) + (count[["C"]] > 1) + (count[["G"]] > 1) + (count[["T"]]  > 1)) > 1)
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
    for (i in 1:numind) {
        if (is.na(ypos[i])) next
        rect(xleft=xmin, xright=xmax,
             ybottom=ypos[i]-0.5, ytop=ypos[i]+0.5, border=NA,
             col=ifelse(x$sites[,i]==majAllele,majColor,ifelse(x$sites[,i]=="N", rgb(0,0,0,0.1), minColor)))
                                        #col=as.character(baseColor[as.character(x$sites[,i])]))
        if (is.null(textColor) || is.null(textColor[[names(x$sites)[i]]])) {
            textcol <- "black"
        } else textcol <- textColor[[names(x$sites)[i]]]
        mtext(names(x$sites)[i], side=2, at=ypos[i], las=2, col=textcol, cex=textSize)
    }
    abline(h=seq(from=min(ypos, na.rm=TRUE)-0.5, to=max(ypos,na.rm=TRUE)+0.5, by=2), lty=3, col="gray")
    invisible(NULL)
}
