write.xls <- function(x, file, col.names=TRUE, ...) write.table(as.data.frame(x), file=file, ..., col.names=col.names, row.names=FALSE, quote=FALSE, sep="\t")
read.dba <- function(dba.file, sel=1:8, fixnames=FALSE, skip=0, ...){ # read DBA restults or any bed format file and return a genomicranges object {{{
    require(GenomicRanges)
    #require(data.table)
    #peaks <- as.data.frame(fread(dba.file))
    require(readr)
    peaks <- as.data.frame(read_tsv(dba.file, skip=skip))
    if(!is.null(sel)) peaks <- peaks[,sel]
    if(fixnames) names(peaks)[1:3] <- c("seqnames", "start", "end")
    idr.peaks <- makeGRangesFromDataFrame(peaks, TRUE, ...)
    return(idr.peaks)
} #}}}
import.narrowpeak <- function(f){ #{{{
    library(rtracklayer)
    import.bed(f, extraCol=c(signalValue='numeric', pValue='numeric', qValue='numeric', peak='integer'))
} #}}}
