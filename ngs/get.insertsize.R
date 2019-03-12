get.insertsize <- function(bam.file, threshold=1000){ #{{{
    library(Rsamtools)
    p1 <- ScanBamParam(what=c("isize"), flag=scanBamFlag(isProperPair=TRUE))
    isize <- scanBam(bam.file, param=p1)[[1]]$isize
    isize <- isize[isize>0 & isize<threshold]
    return(c(mean=mean(isize, na.rm=TRUE), sd=sd(isize, na.rm=TRUE)))
} #}}}
