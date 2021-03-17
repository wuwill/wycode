#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))

get.orgAnn = function(ref.name){  #  {{{
    if(ref.name == 'hg19'){
        return('org.Hs.eg.db')
    }
    if(ref.name %in% c('mm9', 'mm10', 'grcm38')){
        return('org.Mm.eg.db')
    }
    return(FALSE)
} #}}}
get.txdb = function(ref.name){ #  {{{
    if(ref.name == 'hg19'){
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
        return(txdb)
    }
    if(ref.name == 'mm10'){
        library(TxDb.Mmusculus.UCSC.mm10.ensGene)
        txdb = TxDb.Mmusculus.UCSC.mm10.ensGene
        return(txdb)
    }
    if(ref.name == 'mm9'){
        library(TxDb.Mmusculus.UCSC.mm9.knownGene)
        txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
        return(txdb)
      }
    if(ref.name == 'grcm38'){
      txdb = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset='mmusculus_gene_ensembl', taxonomyId=10090)
      return(txdb)
    }
    if(ref.name == 'bdgp5'){
      txdb = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", host='aug2014.archive.ensembl.org', dataset='dmelanogaster_gene_ensembl')
      return(txdb)
    }
    if(ref.name == 'dm6'){
        txdb = makeTxDbFromUCSC(genome='dm6', tablename='ensGene')
        return(txdb)
    }    
    if(ref.name == 'ce10'){
      txdb = makeTxDbFromUCSC(genome='ce10', tablename='ensGene')
      return(txdb)
    }
    return(FALSE)
  } #}}}
fill.gene.symbol = function(peak.df, mygene.df){ #  {{{
  peak.df$symbol = ''
  unique.mygene = mygene.df[!duplicated(mygene.df$query),]
  rownames(unique.mygene) = unique.mygene$query
  for(i in 1:length(peak.df[,1])){
    x = peak.df[i, 'geneId']
    peak.df[i, 'symbol'] = unique.mygene[x, 'symbol']
  }
  return(peak.df)
}#}}}
peak2DF = function (peakfile, header, ...) { #  {{{
    peak.df <- read.delim(peakfile, header = header, comment.char = "#", ...)
    peak.df[, 2] <- peak.df[, 2] + 1
    return(peak.df)
} #}}}
read.peaks = function(fp){ #  {{{
  peaks = readPeakFile(fp, as="GRanges")
  return(peaks)
} #}}}
create.anno <- function(granges, reference, output, txdb=NULL, orgAnn=NULL){ #{{{
    tss <- paste0(output, ".annotated.xls")
    pie <- paste0(output, ".piechart.png")
    tags <- paste0(output, ".avgProfile.png")
    if(is.null(txdb)) { #{{{
        txdb = get.txdb(reference)
        orgAnn = get.orgAnn(reference)
        # check that the reference is valid
        if(typeof(txdb) != 'S4'){
            print(paste(c('Reference genome name not recognized:', reference)))
            quit('no', 1, FALSE)  
        }
    } #}}}

   if(typeof(orgAnn) == 'S4'){
     peak.anno = annotatePeak(granges, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=orgAnn)
   } else {
     peak.anno = annotatePeak(granges, tssRegion=c(-3000, 3000), TxDb=txdb)
   }

   peak.df = as.data.frame(peak.anno)
   gene.symbols = tryCatch(getGenes(unique(peak.df$geneId), fields=c('symbol')), error=function(cond) return(NA))
   if(is.na(gene.symbols) && FALSE){ #{{{
       a <- geneIDs <- unique(peak.df$geneId)
       ngene <- length(geneIDs)
       gene.symbols <- lapply(1:ceiling(ngene/200), function(i) getGenes(a[((i-1)*200+1):min(length(a), 200*i)], fields='symbol')[c('symbol', 'query', '_id')])
       gene.symbols <- do.call(rbind, gene.symbols)
   } #}}}
   if(!is.na(gene.symbols) && 'symbol' %in% colnames(gene.symbols)){
       final.peak.df = fill.gene.symbol(peak.df, gene.symbols)
   }else{
       final.peak.df = peak.df
       print('No gene symbols found for queried geneIds, not adding symbols to annotated peaks file')
   }

   write.table(final.peak.df, file=tss, sep='\t', row.names=F, quote=F)

# Make graphs if output filepaths were passed on the commandline
   if(!is.null(pie)){
       # Create pavis-like pie chart of gene model
       png(pie, height=1350, width=1350, res = 300)
       plotAnnoPie(peak.anno)
       dev.off()
   }

   if(!is.null(tags)){
       # Make read density vs TSS distance graph
       png(tags, height=1350, width=1350)
       p = plotAvgProf2(granges, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
       print(p)
       dev.off()
   }
} #}}}

###  wyang
findMergedBigWig <- function(parent.dir) { #{{{
    bigwigs <- find.file(parent.dir, "*fc.signal.bigwig")
    bigwig.pooled <- grep("call-macs2_pooled/execution", bigwigs, value=TRUE)
    if(length(bigwig.pooled)>=1) return(bigwig.pooled)
    grep("call-macs2/shard-0", bigwigs, value=TRUE)
} #}}}
find.file <- function(path, pattern, pattern2=NULL, vpattern=NULL){ #{{{
    cmd <- paste0("find ", path, " -name '", pattern, "'")
    if(!is.null(pattern2)) cmd <- paste0(cmd, " | grep '", pattern2, "'")
    if(!is.null(vpattern)) cmd <- paste0(cmd, " | grep -v '", vpattern, "'")
    return(system(cmd, TRUE))
} #}}}
chipseq2_annotatePeak <- function(parent.dir,
                         GENOME="hg19", gtf=NULL,
                         peak.files = NULL, bw.files = NULL,
                         chrs=paste0("chr", c(1:30, "X", "Y")),
                         gene_ann=NULL) { #{{{
    library(stringr)
    library(ChIPseeker)
    source("~/wycode/ngs/io.R")
    source("~/wycode/ngs/chipseq.heatmap.R")
    o.dir <- setwd(parent.dir); on.exit(setwd(o.dir))
    qc.files <- find.file(".", "qc.html")
    comparisons <- str_extract(qc.files, '(?<=^./).*(?=/cromwell)')
    if(is.null(peak.files)) peak.files <- sapply(comparisons, function(x) find.file(x, "opti*narrowPeak.gz", "execution/opti"))
    if(is.null(bw.files)) bw.files <- sapply(comparisons, findMergedBigWig)
    groups <- str_replace(comparisons, "_minus_.*$", "")
    if(any(duplicated(groups))) groups <- comparisons
    chip.dir <- data.frame(group=groups, top.path=parent.dir, dir = comparisons, stringsAsFactors = FALSE)
    peaks <- lapply(peak.files, import.narrowpeak)
    names(peaks) <- names(bw.files) <- chip.dir$group

    GENOME <- GENOME
    if(!is.null(gtf)) GTF <- gtf
    source("~/wycode/ngs/chipann.R")
    if(!is.null(gene_ann)) {
        geneAnn <- gene_ann
        TSS_ <- promoters(geneAnn, 0, 1); names(TSS_) <- TSS_$symbol
    }
    reduced.tss <- reduce(TSS_)

    dir.create("peak.annotation")
    #ann.peaks <- lapply(seq_along(chip.dir$group), function(i) my.annotatePeak(peaks[[i]], pdf.file=paste0("peak.annotation/peak.function_annotation.", chip.dir$group[i], ".pdf")))

    chrs <- chrs[chrs %in% seqnames(reduced.tss)]
    reduced.tss <- reduced.tss[seqnames(reduced.tss) %in% chrs]
    seqlevels(reduced.tss) <- chrs
    if(TRUE) { #{{{
    tss.aggregate.signal <- lapply(bw.files, function(x) read.bigwig2rle(bigwig=x, peaks=reduced.tss, bp=1000, aggregate=TRUE))
    names(tss.aggregate.signal) <- names(bw.files)
    save(tss.aggregate.signal, file="peak.annotation/tss.aggregate.signal.RData")

    tss.aggregate.signal <- sapply(tss.aggregate.signal, I)
    library(paletteer)
    library(ggplot2)
    #cls <- paletteer_dynamic("cartography::multi.pal", nrow(chip.dir))
    cls <- c("#CA4D91FF", "#88D651FF", "#8C49CAFF", "#D2B443FF", "#4A3265FF", "#68D1A0FF", "#D05138FF", "#79B0C4FF", "#6A322DFF", "#C8CDA3FF", "#9D8CCBFF", "#617E37FF", "#C38878FF", "#3D4F40FF")[1:length(comparisons)]
    cls <- alpha(cls, 0.6)
    names(cls) <- chip.dir$group

    pdf("peak.annotation/chip_profile.around_TSS.pdf", width=4.5, height=4.5)
    #for(pattern in c("K4me1", "K27", "Mll", "CTCF", "CBP", "300")){
        #groups <- #grep(pattern, chip.dir$group, value=TRUE)
        groups <- names(cls)
        my.mk.signal.histgram(tss.aggregate.signal[, groups], col=cls[groups])
    #}
    dev.off()
    } #}}}


    library(EnrichedHeatmap)
    library(circlize)
    heat.mats <- lapply(bw.files, function(x) get.heat.mat(bigwig = x, peaks=reduced.tss, bp=1000))
    names(heat.mats) <- names(bw.files)
    save(heat.mats, file="peak.annotation/tss.heatmap.matrix.RData")
    strand <- strand(reduced.tss)
    minus <- which(strand %in% "-")
    pdf("peak.annotation/heatmap.around_TSS.minus_reversed.pdf", width=3, height=6)
    for(i in seq_along(heat.mats)){
        mat1 <- heat.mats[[i]]
        mat2 <- mat1[minus,,drop=FALSE]
        mat1[minus,] <- mat2[, ncol(mat2):1]
        name <- names(heat.mats)[i]
        col_fun = colorRamp2(quantile(mat1, c(0, 0.98)), c("white", "red"))
        print(EnrichedHeatmap(mat1, col = col_fun, name = name, ))
    }
    dev.off()

# peak functions
    for(i in seq_along(peaks)){
        ann <- annotatePeakInBatch(peaks[[i]], AnnotationData = geneAnn, multiple = FALSE)
        write_tsv(as.data.frame(ann), paste0("peak.annotation/", names(peaks)[i], ".optimal_set_peak.annotated.xls"))
    }
    for(i in seq_along(peaks)){
        name <- names(peaks)[i]
        create.anno(peaks[[i]], GENOME, paste0("peak.annotation/", name))
    }
    #source("~/wycode/ngs/find.file.R")
    return(list(peaks=peaks, bw.files = bw.files))
} #}}}
