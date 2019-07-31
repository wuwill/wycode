write.xls <- function(x, file, col.names=TRUE, ...) write.table(as.data.frame(x), file=file, ..., col.names=col.names, row.names=FALSE, quote=FALSE, sep="\t")
read.dba <- function(dba.file, sel=NULL, fixnames=FALSE, skip=0, ...){ # read DBA restults or any bed format file and return a genomicranges object {{{
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
read.pos2gr <- function(file, head=TRUE, sel=TRUE, sep="\t", ...){ #{{{
    library(GenomicRanges)
    d <- read.table(file, head=head, sep=sep, ..., as.is=TRUE)[,sel,drop=FALSE]
    d <- data.frame(seqnames=d[[1]], start=d[[2]], end=d[[2]], d[,-(1:2), drop=FALSE], stringsAsFactors=FALSE)
    makeGRangesFromDataFrame(d, keep.extra.columns=TRUE)
} #}}}
read.range2gr <- function(file, head=TRUE, sel=TRUE, sep="\t", ...){ #{{{
    library(GenomicRanges)
    d <- read.table(file, head=head, sep=sep, ..., as.is=TRUE)[,sel,drop=FALSE]
    d <- data.frame(seqnames=d[[1]], start=d[[2]], end=d[[3]], d[,-(1:3), drop=FALSE], stringsAsFactors=FALSE)
    makeGRangesFromDataFrame(d, keep.extra.columns=TRUE)
} #}}}
get.readCounts <- function(bam_file, gr){ #{{{
    library("VariantAnnotation")
    library("Rsamtools")
    sbp = ScanBamParam(which=gr)
    pup = PileupParam(max_depth=1000,
                      min_base_quality=13,
                      min_mapq=0,
                      min_nucleotide_depth=1,
                      min_minor_allele_depth=0,
                      distinguish_strands=FALSE,
                      distinguish_nucleotides=TRUE,
                      ignore_query_Ns=TRUE,
                      include_deletions=FALSE,
                      include_insertions=FALSE)
    x2 <- pileup(bam_file, scanBamParam=sbp, pileupParam=pup)
    x2$allele_pos <- paste0(x2$which_label, "::", x2$nucleotide)
    gr$ref_pos <- paste0(seqnames(gr), ":", start(gr), "-", end(gr), "::", gr$ref)
    gr$alt_pos <- paste0(seqnames(gr), ":", start(gr), "-", end(gr), "::", gr$var)
    gr$n_ref <- gr$n_alt <- 0
    pos <- intersect(gr$ref_pos, x2$allele_pos); gr$n_ref[match(pos, gr$ref_pos)] <- x2$count[match(pos, x2$allele_pos)]
    pos <- intersect(gr$alt_pos, x2$allele_pos); gr$n_alt[match(pos, gr$alt_pos)] <- x2$count[match(pos, x2$allele_pos)]
    gr$vaf <- round((gr$n_alt + 1e-10) / (gr$n_ref + gr$n_alt + 2e-10) * 100, 2)
    return(gr)
} #}}}
