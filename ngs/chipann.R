library(GenomicFeatures)
library(ChIPpeakAnno)
library(ChIPseeker)
#GENOME <- "h19"
REF_DIR <- if(file.exists("/scratch/ref")) "/scratch/ref/gtac/reference_sequences" else "~/local"
if(exists("GENOME")){ #{{{
    if(GENOME == "hg19") {
        library(org.Hs.eg.db)
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        orgAnn <- "org.Hs.eg.db"
        load(file.path(REF_DIR, "chipseq_pipeline_genome_data/hg19.geneAnn.RData"))
    }
    if(GENOME == "mm10"){
        library(org.Mm.eg.db)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        orgAnn <- "org.Mm.eg.db"
        load(file.path(REF_DIR, "chipseq_pipeline_genome_data/mm10.geneAnn.RData"))
    }
    if(GENOME == "mm9"){
        library(org.Mm.eg.db)
        library(TxDb.Mmusculus.UCSC.mm9.knownGene)
        TxDb <- TxDb.Mmusculus.UCSC.mm9.knownGene
        orgAnn <- "org.Mm.eg.db"
        load(file.path(REF_DIR, "chipseq_pipeline_genome_data/mm9.geneAnn.RData"))
    }
} #}}}
read.dba <- function(dba.file, sel=1:8, fixnames=FALSE, skip=0, ...){ # read DBA restults and return a genomicranges object {{{
    #require(data.table)
    require(readr)
    require(GenomicRanges)
    #peaks <- as.data.frame(fread(dba.file))
    peaks <- as.data.frame(read_tsv(dba.file, skip=skip))
    if(!is.null(sel)) peaks <- peaks[,sel]
    if(fixnames) names(peaks)[1:3] <- c("seqnames", "start", "end")
    idr.peaks <- makeGRangesFromDataFrame(peaks, TRUE, ...)
    return(idr.peaks)
} #}}}narrowPeak.hammock.gz
write.xls <- function(x, col.names=TRUE, ...) write.table(as.data.frame(x), ..., col.names=col.names, row.names=FALSE, quote=FALSE, sep="\t")
import.narrowpeak <- function(f){ #{{{
    library(rtracklayer)
    import.bed(f, extraCol=c(signalValue='numeric', pValue='numeric', qValue='numeric', peak='integer'))
} #}}}
getGeneAnn <- function(){ # used to generate geneAnn for new species; hg and mm geneAnn already genereated and saved in /scatch/ref {{{
    annoData <- genes(TxDb)
    annoData$feature <- annoData$gene_id
    annoData <- addGeneIDs(annoData, orgAnn='org.Mm.eg.db', feature_id_type='entrez_id', IDs2Add='symbol')
    return(annoData)
} #}}}
annotateByRegion <- function(){ #{{{
    annoData <- genes(TxDb)
    annoData$feature <- annoData$gene_id
    annoData <- addGeneIDs(annoData, orgAnn=orgAnn, feature_id_type='entrez_id', IDs2Add='symbol')
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
    ret$annRegion <- function(peaks, gene=peaks$feature, cl=NULL, correct.strand=FALSE){ # annotate peaks for overlapping with regions in the gene {{{
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
        peaks.by.gene <- by(annotated.peaks, annotated.peaks$feature, ret$annRegion)
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
subsetByOverlapFraction <- function(x, y, fraction=0.5){ #  {{{
    ov <- intersect(x, y)
    ret <- as.data.frame(findOverlaps(x, ov))
    ret$size.x <- width(x)[ret[[1]]]
    ret$size.ov <- width(ov)[ret[[2]]]
    ret$fraction <- ret$size.ov / ret$size.x
    ikeep <- sort(unique(ret[[1]][ret$fraction>=fraction]))
    return(x[ikeep])
} #}}}

relevel.csAnno <- function(anno){ #{{{
    my.level <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
                  "5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron",
                  "Downstream (<=3kb)", "Distal Intergenic")
    if("csAnno" %in% class(anno)) {
        anno@annoStat$Feature <- factor(as.character(anno@annoStat$Feature), levels=my.level)
        return(anno)
    } else {
        lapply(anno, relevel.csAnno)
    }
} #}}}
my.annotatePeak <- function(peak.files, tssRegion=c(-3000, 3000), pdf.file="", ...) { #{{{
    #https://guangchuangyu.github.io/2014/04/visualization-methods-in-chipseeker/
    peakAnno <- lapply(peak.files, annotatePeak, TxDb=TxDb, annoDb=orgAnn, tssRegion=tssRegion, ...)
    peakAnno <- relevel.csAnno(peakAnno)
    if(!is.null(names(peak.files))) names(peakAnno) <- names(peak.files)
    if(length(peak.files)==1) peakAnno <- peakAnno[[1]]
    if(!is.null(pdf.file) && pdf.file!="") {
        pdf(pdf.file); on.exit(dev.off())
        print(plotAnnoBar(peakAnno))
        print(plotDistToTSS(peakAnno))
        if("GRanges" %in% class(peak.files) || length(peak.files)==1){ #{{{
            print(plotAnnoPie(peakAnno))
            print(vennpie(peakAnno))
            print(upsetplot(peakAnno))
            print(upsetplot(peakAnno, vennpie=TRUE))
        }  #}}}
    }
    return(peakAnno)
} #}}}
peakPathway <- function(peakAnno, tssRegion=c(-1000, 1000), fun=c("enrichPathway", "enrichKEGG"), showCategory=15, title="Pathway Enrichment Analysis"){ #{{{
    ## todo: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
    library(ReactomePA) # rectome pathway
    library(DOSE) # disease ontology
    library(clusterProfiler) #  for Gene Ontology and KEGG enrichment 
    if(is.list(peakAnno)){ #{{{
        genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
        if(is.null(genes[[1]]))
            genes = lapply(peakAnnoList, function(i) seq2gene(i, tssRegion=tssRegion, flankDistance=3000, TxDb=TxDb))
        pathway1 <- compareCluster(geneCluster   = genes,
                                   fun           = fun,
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH")
        dotplot(pathway1, showCategory = showCategory, title = title)
        return(pathway1)
    } #}}}

    # if peakAnno is not a list
    gene <- as.data.frame(peakAnno)$geneID
    if(is.null(gene)) gene <- seq2gene(peakAnno, tssRegion = tssRegion, flankDistance = 3000, TxDb=TxDb)
    pathway1 <- get(fun)(gene)
    dotplot(pathway)
    return(pathway)
} #}}}


