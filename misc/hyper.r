my.hyper <- function(i1, i2=NULL, a1, a2=NULL, i0=NULL, a0=NULL){ #{{{
    if(is.null(i0)) i0 <- i1 + i2
    if(!is.null(a0)) {
        if(is.null(a2)) a2 <- a0 - a1 else {
            a3 <- a1 + a2
            a1 <- a1 * a0 / a3
            a2 <- a2 * a0 / a3
        }
    }
    log10p <- phyper(i1, a1, a2, i0, log.p=TRUE, lower.t=FALSE) / log(10)
    #if(log10p>0) log10p <- phyper(i2, b * a2/(a1+a2), b*a1/(a1+a2), i1+i2, log.p=TRUE, lower.t=TRUE) / log(10) + log(2, 10)
    ret1 <- floor(log10p)
    ret2 <- log10p - ret1
    #print(c(10^ret2, ret1))
    return(sprintf("%.2fe%d", 10^ret2, ret1))
} #}}}
