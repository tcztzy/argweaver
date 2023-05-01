##' Make a color transparent
##'
##' @param col The color (as a character vector)
##' @param alpha The degree of transparency (0=invisible, 1=solid).
##' @return A character vector with a transparent color matching original color in RGB values
##' @export
makeTransparent <- function(col, alpha=0.5) {
    tmp <- col2rgb(col)
    rgb(tmp["red",], tmp["green",], tmp["blue",], alpha=255*alpha, maxColorValue=255)
}


##' Get haplotype names from individual names
##' @param inds character string of diploid individual names
##' @return character vector twice as long as inds, with "_1" and "_2" appended
##' to each individual name
##' @export
getHaps <- function(inds) {
    sort(c(sprintf("%s_1", inds), sprintf("%s_2", inds)))
}




##' Probability of going from "a" lineages to "b" lineages
##' @param tOver2N time (in generations) divided by haploid popsize
##' @param a The starting number of lineages
##' @param b The ending number of lineages (should be <= a)
##' @return Probability under coalescenct
##' @export
probCoalCounts <- function(tOver2N, a, b) {
    C <- prod((b:(b+b-1))*(a:(a-b+1))/(a:(a+b-1)))
    s <- exp(-b*(b-1)*tOver2N)
    C2 <- 1
    C3 <- 1
    if (b < a) {
        for (k in (b+1):a) {
            k1 <- k-1
            C2 <- C2 * (b+k1)*(a-k1)/(a+k1)/(k-b)
            C3 <- -C3
            s <- s + C3*exp(-k*k1*tOver2N)*(2*k-1)/(k1+b)*C2
        }
    }
    s <- s/prod(1:b)
    s*C
}



## same as above but using translation of C code
probCoalCounts2 <- function(tOver2N, a, b) {
    lnC1 <-  0.0
    lnC2 <- 0.0
    C3 <- 1.0;

    for (y in 0:(b-1))
        lnC1 = lnC1 + log((b+y)*(a-y)/(a+y))
    s <-  -b*(b-1)*t/2.0/n;

    for (k in (b+1):a) {
        k1 = (k - 1);
        lnC2 <- lnC2 +  log((b+k1)*(a-k1)/(a+k1)/(k-b))
        C3 <- C3 * -1.0
        val <- -k*k1*t/2.0/n + lnC2 + log((2.0*k-1)/(k1+b))
        if (TRUE) {
            loga = s;
            logc = val;
            if (logc > loga) {
                tmp = logc;
                logc = loga;
                loga = tmp;
            }
            s = loga + log(1.0 + C3*exp(logc - loga));
        }
    }

    for (i in 2:b)
        s = s - log(i);
    return(s + lnC1)
}
