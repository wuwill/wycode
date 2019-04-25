args <- commandArgs(trailingOnly=TRUE)
print(args)
#args <- c("call-macs2_pooled/execution/epic2.results.txt", "call-macs2_ppr1/execution/epic2.results.txt", "call-macs2_ppr2/execution/epic2.results.txt", "call-macs2_pooled/execution/epic2.naive_overlaps.txt1")
subsetByOverlapFraction <- function(x, y, fraction=0.5){ #  {{{
    ov <- intersect(x, y)
    ret <- as.data.frame(findOverlaps(x, ov))
    ret$size.x <- width(x)[ret[[1]]]
    ret$size.ov <- width(ov)[ret[[2]]]
    ret$fraction <- ret$size.ov / ret$size.x
    ikeep <- sort(unique(ret[[1]][ret$fraction>=fraction]))
    return(x[ikeep])
} #}}}
read.dba <- function(dba.file, sel=1:8, fixnames=FALSE, ...){ # read DBA restults and return a genomicranges object {{{
    #require(data.table)
    require(readr)
    require(GenomicRanges)
    peaks <- data.frame(read_tsv(dba.file))
    if(!is.null(sel)) peaks <- peaks[,sel]
    if(fixnames) names(peaks)[1:3] <- c("seqnames", "start", "end")
    idr.peaks <- makeGRangesFromDataFrame(peaks, TRUE, ...)
    return(idr.peaks)
} #}}}narrowPeak.hammock.gz
write.xls <- function(x, col.names=TRUE, ...) write.table(as.data.frame(x), ..., col.names=col.names, row.names=FALSE, quote=FALSE, sep="\t")

if(length(args)<3) {
    stop(paste0("At least 3 parameters are required. \n Rscript %cmd% peaks1 peaks2 output"))
}
f.pooled <- args[1]
f.out <- args[length(args)]
f.filter <- setdiff(args, c(f.pooled, f.out))
pooled <- read.dba(f.pooled, sel=NULL, fixnames=TRUE)
filter <- lapply(f.filter, read.dba, sel=NULL, fixnames=TRUE)
ret <- Reduce(subsetByOverlapFraction, filter, init=pooled)
#system(paste0("head -n 1 ", f.pooled, " > ", f.out))
cat("#", file=f.out)
write.xls(ret, file=f.out, append=TRUE, col.names=TRUE)
