library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
read.dba <- function(dba.file, sel=1:8, ...){ # read DBA restults and return a genomicranges object {{{
    require(data.table)
    require(GenomicRanges)
    peaks <- data.frame(fread(dba.file))
    if(!is.null(sel)) peaks <- peaks[,sel]
    idr.peaks <- makeGRangesFromDataFrame(peaks, TRUE, ...)
    return(idr.peaks)
} #}}}narrowPeak.hammock.gz
import.narrowpeak <- function(f){ #{{{
    library(rtracklayer)
    import.bed(f, extraCol=c(signalValue='numeric', pValue='numeric', qValue='numeric', peak='integer'))
} #}}}
getGeneAnn <- function(){ #{{{
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annoData <- genes(TxDb)
    annoData$feature <- annoData$gene_id
    annoData <- addGeneIDs(annoData, orgAnn='org.Mm.eg.db', feature_id_type='entrez_id', IDs2Add='symbol')
    return(annoData)
} #}}}
annotateByRegion <- function(){ #{{{
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annoData <- genes(TxDb)
    annoData$feature <- annoData$gene_id
    annoData <- addGeneIDs(annoData, orgAnn='org.Mm.eg.db', feature_id_type='entrez_id', IDs2Add='symbol')
    proximal.promoter.cutoff=1000
    immediate.downstream.cutoff=1000
    introns <- intronsByTranscript(TxDb)
    fiveUTRs <- fiveUTRsByTranscript(TxDb)
    threeUTRs <- threeUTRsByTranscript(TxDb)
    exons <- exonsBy(TxDb, by="gene")
    cds <- cdsBy(TxDb, by="gene")
    transcripts <- transcriptsBy(TxDb, by="gene")
    tss <- resize(transcripts, width=1, fix="start")

    options(warn = -1)
    try({
        promoters <- promoters(TxDb, upstream = proximal.promoter.cutoff, 
                               downstream = 0)
        distal.upstream <- promoters(TxDb, upstream = 100000, downstream = 0)
        s <- strand(distal.upstream)
        ipos <- which(s %in% "+")
        ineg <- which(s %in% "-")
        end(distal.upstream)[ipos] <- end(distal.upstream)[ipos] - proximal.promoter.cutoff
        start(distal.upstream)[ineg] <- start(distal.upstream)[ineg] + proximal.promoter.cutoff

        immediateDownstream <- flank(transcripts(TxDb), 
                                     width = immediate.downstream.cutoff, start = FALSE, 
                                     use.names = FALSE)
        distal.downstream <- flank(transcripts(TxDb), 
                                     width = 100000, start = FALSE, 
                                     use.names = FALSE)
        s <- strand(distal.downstream)
        ipos <- which(s %in% "+")
        ineg <- which(s %in% "-")
        start(distal.downstream)[ipos] <- start(distal.downstream)[ipos] + immediate.downstream.cutoff
        end(distal.downstream)[ineg] <- end(distal.downstream)[ineg] - immediate.downstream.cutoff
    })
    my.overlap <- function(x, y) { #{{{
        if(is.null(y)) return(FALSE)
        return(overlapsAny(x, y))
    } #}}}

    ret <- list()
    ret$annRegion <- function(peaks, gene=peaks$feature, cl=NULL, correct.strand=FALSE){ # annoate peaks for overlapping with regions in the gene {{{
        gene <- setdiff(as.character(gene), NA)
        if(length(gene)<1) return(peaks)
        if(length(gene)>1) { #{{{
            ncl <- length(cl)
            n <- length(gene)
            bsize <- 100
            nbatch <- ceiling(n/bsize)
            ret1 <- list()
            iend <- 0
            for(ibatch in 1:nbatch){ #{{{
                if(nbatch>1) cat('\r', iend, "/", n)
                istart <- bsize * (ibatch - 1) + 1
                iend <- min(n, bsize * ibatch)
                reti <- if(ncl>1) parLapply(cl, gene[istart:iend], function(x) ret$annRegion(peaks[peaks$feature %in% x], gene=x)) else
                    lapply(gene[istart:iend], function(x) ret$annRegion(peaks[peaks$feature %in% x], gene=x))
                ret1 <- c(ret1, reti)
            } #}}}
            if(nbatch>1) cat('\n')
            #ret1 <- lapply(gene, function(x) ret$annRegion(peaks[peaks$feature %in% x], gene=x))
            #ret1 <- do.call(getMethod(c, "GenomicRanges"), ret1)
            ret1 <- do.call("c", ret1)
            return(ret1)
        } #}}}

        txs <- as.character(transcripts[[gene]]$tx_id)
        peaks$tss <- my.overlap(peaks, tss[[gene]])
        d.tss <- as.data.frame(distanceToNearest(peaks, tss[[gene]]))
        peaks$distanceToTSS <- d.tss$distance
        if(correct.strand){ #{{{
            strand(peaks) <- strand(tss[[gene]])[d.tss$subjectHits[1]]
        } #}}}
        peaks$threeUTR <- my.overlap(peaks, threeUTRs[names(threeUTRs) %in% txs])
        peaks$fiveUTR <- my.overlap(peaks, fiveUTRs[names(fiveUTRs) %in% txs])
        peaks$exon <- my.overlap(peaks, exons[gene])
        peaks$intron <- my.overlap(peaks, introns[names(introns) %in% txs])
        peaks$promoter <- my.overlap(peaks, promoters[promoters$tx_id %in% txs])
        peaks$distal.upstream <- my.overlap(peaks, distal.upstream[distal.upstream$tx_id %in% txs])
        peaks$immediateDownstream <- my.overlap(peaks, immediateDownstream[immediateDownstream$tx_id %in% txs])
        peaks$distal.downstream <- my.overlap(peaks, distal.downstream[distal.downstream$tx_id %in% txs])
        return(peaks)
    } #}}}
    ret$annByGene <- function(genes, annotated.peaks, use.id=FALSE){ #{{{
        # what peaks overlap the gene?
        annotated.peaks <- if(use.id) annotated.peaks[annotated.peaks$feature %in% genes, ] else
            annotated.peaks[annotated.peaks$symbol %in% genes, ]
        peaks.by.gene <- by(annotated.peaks, annoated.peaks$feature, ret$annRegion)
        return(peaks.by.genen)
    } #}}}
    ret$annPeak <- function(idr.peaks, maxgap=100000){ # annotate IDR peaks to genes with 10kbp {{{
        peak.annotated = annotatePeakInBatch(idr.peaks, AnnotationData=annoData, output='overlapping', maxgap=maxgap)
        peak.annotated <- addGeneIDs(peak.annotated, orgAnn='org.Mm.eg.db', feature_id_type='entrez_id', IDs2Add='symbol')
# change strand with respect to gene orientation
        peak.annotated$feature_strand[is.na(peak.annotated$feature_strand)] <- "*"
        strand(peak.annotated) <- peak.annotated$feature_strand
        return(peak.annotated)
    } #}}}
    ret$tss <- function(gene=NULL){ #{{{
        if(is.null(gene)) tss else tss[[gene]]
    } #}}}
    ret$exons <- function(gene=NULL){ #{{{
        if(is.null(gene)) exons else exons[[gene]]
    } #}}}
    ret$threeUTRs <- function(gene=NULL){ #{{{
        if(is.null(gene)) return(threeUTRs)
        txs <- as.character(transcripts[[gene]]$tx_id)
        threeUTRs[names(threeUTRs) %in% txs]
    } #}}}
    ret$fiveUTRs <- function(gene=NULL){ #{{{
        if(is.null(gene)) return(fiveUTRs)
        txs <- as.character(transcripts[[gene]]$tx_id)
        fiveUTRs[names(fiveUTRs) %in% txs]
    } #}}}
    ret$introns <- function(gene=NULL){ #{{{
        if(is.null(gene)) return(introns)
        txs <- as.character(transcripts[[gene]]$tx_id)
        introns[names(introns) %in% txs]
    } #}}}
    ret$promoters <- function(gene=NULL){ #{{{
        if(is.null(gene)) return(promoters)
        txs <- as.character(transcripts[[gene]]$tx_id)
        promoters[promoters$tx_id %in% txs]
    } #}}}
    ret$immediateDownstream <- function(gene=NULL){ #{{{
        if(is.null(gene)) return(immediateDownstream)
        txs <- as.character(transcripts[[gene]]$tx_id)
        immediateDownstream[immediateDownstream$tx_id %in% txs]
    } #}}}
    ret$get.annoData <- function() return(annoData)

    return(ret)
} #}}}
