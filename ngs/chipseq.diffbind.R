### old way
# need module load bedtools
mk.compare <- function(itest, idr=FALSE, intersect=FALSE){ #{{{
    test.nm <- names(test.groups)[itest]
    test.grs <- test.groups[[itest]]
    if(idr) test.nm <- paste0(test.nm, ".idr")
    if(idr && intersect) test.nm <- paste0(test.nm, ".intersect")
    if(!file.exists(test.nm)) dir.create(test.nm)
    o.dir <- setwd(test.nm); on.exit(setwd(o.dir))
    cat(paste(test.grs, collapse="\t"), "\n", file="compare.txt", sep="")
} #}}}
mk.dba_sample <- function(itest, idr=FALSE, intersect=FALSE){ #{{{
    test.nm <- names(test.groups)[itest]
    test.grs <- test.groups[[itest]]
    if(idr) test.nm <- paste0(test.nm, ".idr")
    if(idr && intersect) test.nm <- paste0(test.nm, ".intersect")
    if(!file.exists(test.nm)) dir.create(test.nm)
    o.dir <- setwd(test.nm); on.exit(setwd(o.dir))
    out <- c("dba_sample.csv")
    get.info <- function(group){ #{{{
        nm <- sample.info$V3[sample.info$V1 %in% group]
        atac.dir <- paste0("../", group, "_atac/out")
        bam.pattern <- paste0(".R1.trim.PE2SE.nodup.bam")
        bams <- system(paste0("ls ", file.path(atac.dir, "align/rep*/*"), bam.pattern), TRUE)
        nm2 <- gsub(bam.pattern, "", basename(bams))
        peak.pattern <- paste0(".R1.trim.PE2SE.nodup.tn5.pf.filt.narrowPeak.gz")
        idr.peak <- system(paste0("ls ", file.path(atac.dir, "peak/idr/optimal*/*narrowPeak.gz")), TRUE)
        peaks <- system(paste0("ls ", file.path(atac.dir, "peak/macs2/rep*/*"), peak.pattern), TRUE)
        nm3 <- gsub(peak.pattern, "", basename(peaks))
        peaks <- peaks[match(nm2, nm3)]
        if(idr) {
            if(intersect){
                # create intersct files
                peaks0 <- peaks
                peaks <- gsub(".gz$", "", basename(peaks0))
                for(i in seq_along(peaks0)){
                    system(paste0("intersectBed -a ", peaks0[i], " -b ", idr.peak, " > ", peaks[i]))
                }
            } else peaks <- idr.peak
            #peaks <- system(paste0("ls ", file.path(atac.dir, "peak/idr/optimal*/*narrowPeak.gz")), TRUE)
        }
        ret <- data.frame(SampleID=nm2, Condition=group, Replicate=seq_along(nm2), bamReads=bams, Peaks=peaks, PeakCaller="narrow", stringsAsFactors=FALSE)
        return(ret)
    } #}}}
    write.table(rbind(get.info(test.grs[1]), get.info(test.grs[2])),
        file="dba_samples.csv", sep=",",
        col.names=TRUE, row.names=FALSE, quote=FALSE)
} #}}}
mk.dba_sample.encode_atac <- function(itest, idr=FALSE, intersect=FALSE){ #{{{
    test.nm <- names(test.groups)[itest]
    test.grs <- test.groups[[itest]]
    if(idr) test.nm <- paste0(test.nm, ".idr")
    if(idr && intersect) test.nm <- paste0(test.nm, ".intersect")
    if(!file.exists(test.nm)) dir.create(test.nm)
    o.dir <- setwd(test.nm); on.exit(setwd(o.dir))
    out <- c("dba_sample.csv")
    get.info <- function(group){ #{{{
        nm <- sample.info$V3[sample.info$V1 %in% group]
        bam.pattern <- paste0("_R1.trim.nodup.bam")
        i.atac <- match(group, atac.dir$group)
        atac.top.dir <- file.path(atac.dir$top.path[i.atac],  atac.dir$dir[i.atac])
        bams <- find.file(atac.top.dir,
                          pattern=paste0("*", bam.pattern),
                          pattern2="call-filter/shard")
        nm2 <- gsub(bam.pattern, "", basename(bams))

        peak.pattern <- "_R1.trim.nodup.*bfilt.narrowPeak.gz"
        idr.peak <- find.file(atac.top.dir, "optimal*narrowPeak.gz", "idr/execution")
        peaks <- find.file(atac.top.dir, paste0("*", peak.pattern), "macs2/shard")
        nm3 <- gsub(peak.pattern, "", basename(peaks))
        peaks <- peaks[match(nm2, nm3)]
        if(idr) {
            if(intersect){
                # create intersct files
                peaks0 <- peaks
                peaks <- gsub(".gz$", "", basename(peaks0))
                for(i in seq_along(peaks0)){
                    system(paste0("intersectBed -a ", peaks0[i], " -b ", idr.peak, " > ", peaks[i]))
                }
            } else peaks <- idr.peak
            #peaks <- system(paste0("ls ", file.path(atac.dir, "peak/idr/optimal*/*narrowPeak.gz")), TRUE)
        }
        ret <- data.frame(SampleID=nm2, Condition=group, Replicate=seq_along(nm2), bamReads=bams, Peaks=peaks, PeakCaller="narrow", stringsAsFactors=FALSE)
        return(ret)
    } #}}}
    write.table(rbind(get.info(test.grs[1]), get.info(test.grs[2])),
        file="dba_samples.csv", sep=",",
        col.names=TRUE, row.names=FALSE, quote=FALSE)
} #}}}
mk.diffbind.job <- function(itest, idr=FALSE, intersect=FALSE){ #{{{
    encode_atac <- exists("atac.dir") && length(find.file(file.path(atac.dir$top.path[1], atac.dir$dir[1]), "align", "out/align")) < 1
    mk.compare(itest, idr, intersect)
    if(encode_atac)
        mk.dba_sample.encode_atac(itest, idr, intersect) else
            mk.dba_sample(itest, idr, intersect)
    test.nm <- names(test.groups)[itest]
    if(idr) test.nm <- paste0(test.nm, ".idr")
    if(idr && intersect) test.nm <- paste0(test.nm, ".intersect")
    o.dir <- setwd(test.nm); on.exit(setwd(o.dir))
    mv <- function(idr=FALSE) { #{{{
        vs <- paste(test.groups[[itest]], collapse="_vs_")
        #input <- file.path("diffbind", paste(c("DBA_", "merged_", "DBA_"), vs, c(".csv", "_peaks.png", "tss.txt"), sep=""))
        input <- paste(c("DBA_", "merged_", "DBA_"), vs, c(".csv", "_peaks.png", "tss.txt"), sep="")
        out <- paste(c("DBA_", "merged_", "DBA_"), vs, c(".idr.csv", "_peaks.idr.png", "annotated.idr.xls"), sep="")
        ret <-paste0("mv ", input,  " ", out, "\n")
        if(idr) return(paste(ret, collapse="")) else return(ret[3])
    } #}}}
    cat("#!/bin/bash -e
#SBATCH --export=ALL
#SBATCH --mem=12000
cd ", getwd(), "

Rscript /scratch/gtac/software/atac/diff_bind.R --reference ", genome, "\n",
mv(idr), file="runDiffBind.sh")
    system("sbatch --mem=12000 runDiffBind.sh")
#Rscript /scratch/gtac/software/python_virtual-env/chipseq/bin/diff_bind.R --reference ", genome, "\n",
#Rscript /scratch/gtac/software/apps/chipseq/scripts/diff_bind.R --reference ", genome, "\n",
#Rscript /scratch/gtac/software/apps/miniconda3/envs/chipseq/bin/diff_bind.R --reference ", genome, "\n",
} #}}}
#sample.info <- data.frame(V1=group, V2=index, V3=name, stringsAsFactors=FALSE)
#test.groups <- list(HPC_WTvsHSC_WT=c("HPC_WT", "HSC_WT"))
#sapply(seq_along(test.groups), mk.diffbind.job, idr=TRUE, intersect=TRUE)

### new way testing any regions
# source find.file.R / io.R
get.atac.sampleSheet <- function(group, atac.dir, ...){ #{{{
    ret <- list()
    i <- match(group, atac.dir$group)
    path <- file.path(atac.dir$top.path[i], atac.dir$dir[i], "out")
    samples <- list.files(file.path(path, "align"), "^rep")
    tagAligns <- sapply(file.path(path, "align", samples, "*tn5.tagAlign.gz"), Sys.glob)
    beds <- gsub("tagAlign.gz", "bed.gz", tagAligns)
    file.symlink(tagAligns, beds)
    peaks <- sapply(file.path(path, "peak/macs2", samples, "*filt.narrowPeak.gz"), Sys.glob)
    sampleID <- paste0(group, "_", samples)
    ctrl <- file.path(path, "align/ctl")
    ctrl <- if(file.exists(ctrl)) "ctrl" else ""
    if(ctrl=="") {
        bedctl <- ""
    } else {
        ctrl.tagAligns <- Sys.glob(file.path(path, "align", ctrl, "*tn5.tagAlign.gz"))
        bedctl <- gsub("tagAlign.gz", "bed.gz", ctrl.tagAligns)
        file.symlink(ctrl.tagAligns, bedctl)
    }
    ret <- data.frame(SampleID=sampleID, Tissue="", Factor=group, Condition=group, Replicate=gsub("rep", "", samples),
                      bamReads=beds, bamControl=bedctl, Peaks=peaks, PeakCaller="macs", stringsAsFactors=FALSE)
    return(ret)
} #}}}
get.chip2.sampleSheet <- function(group, chip.dir, peak.file="optimal", useBam=FALSE){ #{{{
    #peak.file: NULL - find for each sample; "optimal" optimal peaks; or self specified peaks
    ret <- list()
    i <- match(group, chip.dir$group)
    path <- file.path(chip.dir$top.path[i], chip.dir$dir[i])
    #samples <- list.files(file.path(path, "align"), "^rep")
    #tagAligns <- sapply(file.path(path, "align", samples, "*tn5.tagAlign.gz"), Sys.glob)
    if(useBam){ #{{{
        tagAligns <- beds <- find.file(path, "*.nodup.bam", pattern2="filter/")
        samples <- sapply(strsplit(tagAligns, "/"), grep, pattern="^shard", value=TRUE)
    } else {
        tagAligns <- find.file(path, "*nodup.*tagAlign.gz", pattern2="bam2ta/")
        samples <- sapply(strsplit(tagAligns, "/"), grep, pattern="^shard", value=TRUE)
        beds <- gsub("tagAlign.gz", "bed.gz", tagAligns)
        file.symlink(tagAligns, beds)
    }#}}}
    if(is.null(peak.file)){ #{{{
        peaks <- find.file(path, "*bfilt.narrowPeak.gz", pattern2="macs2/")
    } else if(peak.file==c("optimal", "idr", "overlap")) {
        peaks <- find.file(path, "*optimal*narrowPeak.gz", pattern2="/execution") 
        if(length(peaks)>1) peaks <- grep(peak.file, peaks, value=TRUE)
        #if(length(peaks)>1) peaks <- find.file(path, "*optimal*narrowPeak.gz", pattern2="idr")
    } else {
        peaks <- peak.file
    } #}}}
    sampleID <- paste0(group, "_", samples)
    if(useBam) {  #  {{{
        pooled.bam <- find.file(path, "pooled.merged.nodup.bam", pattern2="filter_ctl")
        all.bams <- find.file(path, "*nodup.bam", pattern2="filter_ctl")
        n.bams <- length(all.bams)
        if(length(pooled.bam) == 0){
            if(n.bams==0) {
                bedctl <- ctrl <- ""
            } else if(n.bams==1) {
                bedctl <- all.bams
            } else {
                library(Rsamtools)
                destBam <- gsub("call-filter_ctl/.*$", "call-filter_ctl/pooled.merged.nodup.bam", all.bams[1])
                mergeBam(all.bams, destBam)
                indexBam(destBam)
                bedctl <- destBam
            }
        } else {
            bedctl <- pooled.bam
        } #}}}
    }  else { #  {{{
        ctrl.tagAligns <- find.file(path, "*nodup.pooled.tagAlign.gz", pattern2="pool_ta_ctl")
        if(length(ctrl.tagAligns)==0) {
            bedctl <- ctrl <- ""
        } else {
            ctrl <- paste0("ctrl.", sapply(strsplit(ctrl.tagAligns, "/"), grep, pattern="^shard", value=TRUE))
            bedctl <- gsub("tagAlign.gz", "bed.gz", ctrl.tagAligns)
            file.symlink(ctrl.tagAligns, bedctl)
        }
    } #}}}
    ret <- data.frame(SampleID=sampleID, Tissue="", Factor=group, Condition=group, Replicate=gsub("shard", "", samples),
                      bamReads=beds, bamControl=bedctl, Peaks=peaks, PeakCaller="macs", stringsAsFactors=FALSE)
    return(ret)
} #}}}
my.grange.intersect <- function(x, y) { # interect of x and y while keeping meta data of x{{{
    ret <- intersect(x, y)
    ov <- as.data.frame(findOverlaps(ret, x))
    mcols(ret) <- mcols(x)[ov[match(1:length(ret), ov[[1]]),2],,drop=FALSE]
    return(ret)
} #}}}
get.sampleSheet <- function(group, chip.dir, peak.file="optimal", useBam=TRUE, diffbind.workdir=NULL, atac=FALSE, intersect.peak=NULL){ #{{{
    #peak.file: NULL - find for each sample; "optimal"/"idr" optimal peaks; "idr.intersect" - peaks for each sample intersect with IDR ; or self specified peaks
    ret <- list()
    i <- match(group, chip.dir$group)
    path <- file.path(chip.dir$top.path[i], chip.dir$dir[i])
    #samples <- list.files(file.path(path, "align"), "^rep")
    #tagAligns <- sapply(file.path(path, "align", samples, "*tn5.tagAlign.gz"), Sys.glob)
    grep.peak.file <- function(..., exact=FALSE) { #{{{
        patterns <- unlist(list(...))
        if(length(peak.file) != 1) return(FALSE)
        if(exact) peak.file %in% c(patterns, paste0(patterns, ".intersect")) else any(sapply(patterns, grepl, x=peak.file))
    } #}}}
    if(is.null(peak.file)){ #{{{
        peaks <- find.file(path, "*bfilt.narrowPeak.gz", pattern2="macs2/")
    } else if(grep.peak.file("optimal", "idr", "overlap", exact=TRUE) || !is.null(intersect.peak)) {
        idr <- if(is.null(intersect.peak)) find.file(path, "*optimal*narrowPeak.gz", pattern2="/execution")  else intersect.peak
        if(length(idr)>1) idr <- grep(gsub(".intersect", "", peak.file), idr, value=TRUE)
        if(grep.peak.file("intersect$")) {
            peaks0 <- find.file(path, "*bfilt.narrowPeak.gz", pattern2="macs2/")
            peaks <- file.path(diffbind.workdir, gsub(".narrowPeak.gz", paste0(".", peak.file, ".narrowPeak.gz"), basename(peaks0)))
            peak.idr <- import.narrowpeak(idr)
            for(i in seq_along(peaks)) {
                peaki <- import.narrowpeak(peaks0[i])
                export.bed(my.grange.intersect(peaki, peak.idr), peaks[i])
            }
        } else {
            peaks <- idr
        }
    } else {
        peaks <- peak.file
    } #}}}
    if(useBam){ #{{{
        tagAligns <- beds <- find.file(path, "*.nodup.bam", pattern2="filter/")
        samples <- sapply(strsplit(tagAligns, "/"), grep, pattern="^shard", value=TRUE)
    } else {
        tagAligns <- find.file(path, "*nodup.*tagAlign.gz", pattern2="bam2ta/")
        samples <- sapply(strsplit(tagAligns, "/"), grep, pattern="^shard", value=TRUE)
        beds <- gsub("tagAlign.gz", "bed.gz", tagAligns)
        file.symlink(tagAligns, beds)
    }#}}}
    sampleID <- paste0(group, "_", samples)
    if(atac) {
        ret <- data.frame(SampleID=sampleID, Tissue="", Factor=group, Condition=group, Replicate=gsub("shard", "", samples),
                          bamReads=beds, Peaks=peaks, PeakCaller="macs", stringsAsFactors=FALSE)

    } else {
        if(useBam) {  #  {{{
            pooled.bam <- find.file(path, "pooled.merged.nodup.bam", pattern2="filter_ctl")
            all.bams <- find.file(path, "*nodup.bam", pattern2="filter_ctl")
            n.bams <- length(all.bams)
            if(length(pooled.bam) == 0){
                if(n.bams==0) {
                    bedctl <- ctrl <- ""
                } else if(n.bams==1) {
                    bedctl <- all.bams
                } else {
                    library(Rsamtools)
                    destBam <- gsub("call-filter_ctl/.*$", "call-filter_ctl/pooled.merged.nodup.bam", all.bams[1])
                    mergeBam(all.bams, destBam)
                    indexBam(destBam)
                    bedctl <- destBam
                }
            } else {
                bedctl <- pooled.bam
            } #}}}
        }  else { #  {{{
            ctrl.tagAligns <- find.file(path, "*nodup.pooled.tagAlign.gz", pattern2="pool_ta_ctl")
            if(length(ctrl.tagAligns)==0) {
                bedctl <- ctrl <- ""
            } else {
                ctrl <- paste0("ctrl.", sapply(strsplit(ctrl.tagAligns, "/"), grep, pattern="^shard", value=TRUE))
                bedctl <- gsub("tagAlign.gz", "bed.gz", ctrl.tagAligns)
                file.symlink(ctrl.tagAligns, bedctl)
            }
        } #}}}
        ret <- data.frame(SampleID=sampleID, Tissue="", Factor=group, Condition=group, Replicate=gsub("shard", "", samples),
                          bamReads=beds, bamControl=bedctl, Peaks=peaks, PeakCaller="macs", stringsAsFactors=FALSE)
    }
    return(ret)
} #}}}
library(DiffBind)
run.DiffBind <- function(group1, group2, atac.dir, peak.file=NULL, atac=TRUE, dba=NULL, useBam=FALSE, out.dir=NULL, sampleSheet=NULL,
                         method=DBA_EDGER, bTagwise=FALSE, useCtl=TRUE, bUseSummarizeOverlaps=FALSE, ...){ #{{{
    n.clean <- 0
    if(is.null(out.dir)) out.dir <- paste0("DBA.", group1, "-", group2)
    if(!file.exists(out.dir)) dir.create(out.dir)
    o.dir <- setwd(out.dir); n.clean <- 1
    on.exit({
                if(n.clean>0) setwd(o.dir)
                if(n.clean>1) dev.off()
            })
    group1.dir <- atac.dir[grep(group1, atac.dir$group),, drop=FALSE]
    if(is.null(dba)){ #{{{
        if(is.null(sampleSheet)) {
            get.sampleSheet <- if(atac){ #{{{
                if(group1.dir$top.path %in% new.pipe.list) get.sampleSheet else get.atac.sampleSheet
            } else {
                if(group1.dir$top.path %in% new.pipe.list) get.sampleSheet else get.chip.sampleSheet
            } #}}}
            sampleSheet <- rbind(get.sampleSheet(group1, atac.dir, peak.file=peak.file, useBam=useBam, atac=atac, diffbind.workdir=getwd()),
                                          get.sampleSheet(group2, atac.dir, peak.file=peak.file, useBam=useBam, atac=atac, diffbind.workdir=getwd()))
            rownames(sampleSheet) <- NULL
            if(min(table(sampleSheet$Condition))<2) {
                cat("Warning: only 1 sample in some groups!\n")
                print(sampleSheet)
                return()
            }
            write.xls(sampleSheet, file="sampleSheet.txt")
        } else if(length(sampleSheet)==1){
            library(data.table)
            sampleSheet <- data.frame(fread(sampleSheet, ...))
        }
        if(!useCtl) {
            cat("Set bamControl to NULL!\n")
            sampleSheet$bamControl <- NULL
        }
        dba <- dba(sampleSheet=sampleSheet)
        #dba <- dba.count(dba, summits=250, bUseSummarizeOverlaps=TRUE)
        dba <- if(useBam) dba.count(dba, bUseSummarizeOverlaps=bUseSummarizeOverlaps) else dba.count(dba, summits=250) 
        dba <- dba.contrast(dba, categories=DBA_CONDITION, minMembers=2)
        dba <- dba.analyze(dba, method=method, bTagwise=bTagwise, bSubControl=TRUE, ...)
		 #analysis = dba.analyze(contrast, method=METHOD, bCorPlot=F, bTagwise=F, bSubControl=T)
        save(dba, file="dba.ret.RData")
    }#}}}
    my.dba.report <- function(...){ #{{{
        ret <- dba.report(..., bCalled=FALSE, method=method)
        if(length(ret)<2) return(ret)
        dba.report(..., bCalled=TRUE, method=method)
    } #}}}
    dba.report1 <- my.dba.report(dba, th=0.1)
    dba.report2 <- my.dba.report(dba, th=0.01, bUsePval=TRUE, fold=2)
    dba.report3 <- my.dba.report(dba, th=0.05)
    out.dir <- basename(out.dir)
    file1 <- paste0("DBA.", out.dir, '.', method, ".FDR0.1.report.xls")
    file2 <- paste0("DBA.", out.dir, '.', method, ".p0.01.Fold2.report.xls")
    file3 <- paste0("DBA.", out.dir, '.', method, ".FDR0.05.report.xls")
    file.pdf <- paste0("DBA.", out.dir, ".pdf")
    write.xls(dba.report1, file=file1)
    write.xls(dba.report2, file=file2)
    write.xls(dba.report3, file=file3)
    pdf(file.pdf); n.clean <- 2
    dba.plotHeatmap(dba, method=method); #dba.plotVenn(dba)
    dba.plotMA(dba, method=method, th=0.1);
    dba.plotMA(dba, bUsePval=TRUE, fold=2, method=method, th=0.01);
    dba.plotMA(dba, method=method, th=0.05);
    if(length(dba.report1)>0) dba.plotVolcano(dba, method=method, th=0.1)
    if(length(dba.report2)>0) dba.plotVolcano(dba, bUsePval=TRUE, fold=2, method=method, th=0.01);
    if(length(dba.report3)>0) dba.plotVolcano(dba, method=method, th=0.05)
    save(dba, dba.report1, dba.report2, dba.report3, file="DBA.RData")
    setwd(o.dir)
} #}}}

#atac.dir <- data.frame(
                       #group=c("BJFF", "KMT2CHET", "KMT2CKO"),
                       #top.path=c(atac.dir0, atac.dir0, atac.dir0),
                       #dir=c("BJFF_atac", "KMT2Chet_atac", "KMT2C-KO_atac")
                       #)
#run.DiffBind("BJFF", "KMT2CKO", atac.dir)
