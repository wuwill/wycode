#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))

get.orgAnn = function(ref.name){
    if(ref.name == 'hg19'){
        return('org.Hs.eg.db')
    }
    if(ref.name %in% c('mm9', 'mm10', 'grcm38')){
        return('org.Mm.eg.db')
    }
    return(FALSE)
}

get.txdb = function(ref.name){
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
  }


fill.gene.symbol = function(peak.df, mygene.df){
  peak.df$symbol = ''
  unique.mygene = mygene.df[!duplicated(mygene.df$query),]
  rownames(unique.mygene) = unique.mygene$query
  for(i in 1:length(peak.df[,1])){
    x = peak.df[i, 'geneId']
    peak.df[i, 'symbol'] = unique.mygene[x, 'symbol']
  }
  return(peak.df)
}

peak2DF = function (peakfile, header, ...)
{
    peak.df <- read.delim(peakfile, header = header, comment.char = "#", ...)
    peak.df[, 2] <- peak.df[, 2] + 1
    return(peak.df)
}

read.peaks = function(fp){
  peaks = readPeakFile(fp, as="GRanges")
  return(peaks)
}

parser = ArgumentParser()
parser$add_argument('--peak', required=T)
parser$add_argument('--tss', required=T)
parser$add_argument('--reference', required=T)
parser$add_argument('--pie')
parser$add_argument('--tags')
parser$add_argument('--sep', default='', help='sep character to parse peak file with')
parser$add_argument('--noheader', action='store_true', help='peaks file has no header row')


# load phony args for interactive debugging
test.args = str_split('--tss test.xls --reference mm10 --peak MACS_peaks_mm10.xls --pie pie.png --tags tag.png', ' ')[[1]]

if(!interactive()) args = parser$parse_args() else args = parser$parse_args(test.args)

txdb = get.txdb(args$reference)
orgAnn = get.orgAnn(args$reference)

print("@@0.3")
# check that the reference is valid
if(typeof(txdb) != 'S4'){
     print(paste(c('Reference genome name not recognized:', args$reference)))
     quit('no', 1, FALSE)  
}

tmp.peak.df = peak2DF(args$peak, !args$noheader, sep=args$sep)
granges = ChIPseeker:::peakDF2GRanges(tmp.peak.df)

# Create annotated peaks
   if(typeof(orgAnn) == 'S4'){
     peak.anno = annotatePeak(granges, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=orgAnn)
   }else{
     peak.anno = annotatePeak(granges, tssRegion=c(-3000, 3000), TxDb=txdb)
   }
print("@@1")

peak.df = as.data.frame(peak.anno)
#gene.symbols = getGenes(unique(peak.df$geneId), fields=c('symbol'))
gene.symbols = tryCatch(getGenes(unique(peak.df$geneId), fields=c('symbol')), error=function(cond) return(NA))
print("@@2")
if(is.na(gene.symbols)){ #{{{
    a <- geneIDs <- unique(peak.df$geneId)
    ngene <- length(geneIDs)
    gene.symbols <- lapply(1:ceiling(ngene/200), function(i) getGenes(a[((i-1)*200+1):min(length(a), 200*i)], fields='symbol')[c('symbol', 'query', '_id')])
    gene.symbols <- do.call(rbind, gene.symbols)
} #}}}
print("@@3")
if('symbol' %in% colnames(gene.symbols)){
    final.peak.df = fill.gene.symbol(peak.df, gene.symbols)
}else{
    final.peak.df = peak.df
    print('No gene symbols found for queried geneIds, not adding symbols to annotated peaks file')
}

print("@@4")
write.table(final.peak.df, file=args$tss, sep='\t', row.names=F, quote=F)

# Make graphs if output filepaths were passed on the commandline
if(!is.null(args$pie)){
  # Create pavis-like pie chart of gene model
  png(args$pie, height=480, width=480)
  plotAnnoPie(peak.anno)
  dev.off()
}

if(!is.null(args$tags)){
  # Make read density vs TSS distance graph
  png(args$tags, height=480, width=480)
  p = plotAvgProf2(granges, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  print(p)
  dev.off()
}
