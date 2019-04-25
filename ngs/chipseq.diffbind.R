### old way
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
mk.diffbind.job <- function(itest, idr=FALSE, intersect=FALSE){ #{{{
    mk.compare(itest, idr, intersect)
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
cd ", getwd(), "

Rscript /scratch/gtac/software/apps/atac/diff_bind.R --reference ", genome, "\n",
mv(idr), file="runDiffBind.sh")
    system("sbatch --mem=12000 runDiffBind.sh")
#Rscript /scratch/gtac/software/python_virtual-env/chipseq/bin/diff_bind.R --reference ", genome, "\n",
#Rscript /scratch/gtac/software/apps/chipseq/scripts/diff_bind.R --reference ", genome, "\n",
#Rscript /scratch/gtac/software/apps/miniconda3/envs/chipseq/bin/diff_bind.R --reference ", genome, "\n",
} #}}}
#test.groups <- list(HPC_WTvsHSC_WT=c("HPC_WT", "HSC_WT"))
#sapply(seq_along(test.groups), mk.diffbind.job, idr=TRUE, intersect=TRUE)

### new way testing any regions
find.file <- function(path, pattern, pattern2=NULL, vpattern=NULL){ #{{{
    cmd <- paste0("find ", path, " -name '", pattern, "'")
    if(!is.null(pattern2)) cmd <- paste0(cmd, " | grep '", pattern2, "'")
    if(!is.null(vpattern)) cmd <- paste0(cmd, " | grep -v '", vpattern, "'")
    return(system(cmd, TRUE))
} #}}}
write.xls <- function(x, ...) write.table(as.data.frame(x), ..., col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
get.atac.sampleSheet <- function(group, atac.dir){ #{{{
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
get.chip2.sampleSheet <- function(group, chip.dir, peak.file="optimal"){ #{{{
    #peak.file: NULL - find for each sample; "optimal" optimal peaks; or self specified peaks
    ret <- list()
    i <- match(group, chip.dir$group)
    path <- file.path(chip.dir$top.path[i], chip.dir$dir[i])
    #samples <- list.files(file.path(path, "align"), "^rep")
    #tagAligns <- sapply(file.path(path, "align", samples, "*tn5.tagAlign.gz"), Sys.glob)
    tagAligns <- find.file(path, "*nodup.tagAlign.gz$", pattern2="bam2ta/")
    samples <- sapply(strsplit(tagAligns, "/"), grep, pattern="^shard", value=TRUE)
    beds <- gsub("tagAlign.gz", "bed.gz", tagAligns)
    file.symlink(tagAligns, beds)
    if(is.null(peak.file)){ #{{{
        peaks <- find.file(path, "*bfilt.narrowPeak.gz$", pattern2="macs2/")
    } else if(peak.file=="optimal") {
        peaks <- find.file(path, "*narrowPeak.gz$", pattern2="optimal") 
    } else {
        peaks <- peak.file
    } #}}}
    sampleID <- paste0(group, "_", samples)
    ctrl.tagAligns <- find.file(path, "*nodup.tagAlign.gz$", pattern2="bam2ta_ctl/")
    if(length(ctrl.tagAligns)==0) {
        bedctl <- ""
        ctrl <- ""
    } else {
        ctrl <- paste0("ctrl.", sapply(strsplit(ctrl.tagAligns, "/"), grep, pattern="^shard", value=TRUE))
        bedctl <- gsub("tagAlign.gz", "bed.gz", ctrl.tagAligns)
        file.symlink(ctrl.tagAligns, bedctl)
    }
    ret <- data.frame(SampleID=sampleID, Tissue="", Factor=group, Condition=group, Replicate=gsub("rep", "", samples),
                      bamReads=beds, bamControl=bedctl, Peaks=peaks, PeakCaller="macs", stringsAsFactors=FALSE)
    return(ret)
} #}}}
run.DiffBind <- function(group1, group2, atac.dir, peak.file=NULL, atac=TRUE){ #{{{
    sampleSheet <- if(atac) rbind(get.atac.sampleSheet(group1, atac.dir), get.atac.sampleSheet(group2, atac.dir)) else
        rbind(get.chip2.sampleSheet(group1, atac.dir, peak.file=peak.file), get.chip2.sampleSheet(group2, atac.dir, peak.file=peak.file))
    rownames(sampleSheet) <- NULL
    if(min(table(sampleSheet$Condition))<2) {
        cat("XXX: only 1 sample in some groups!\n")
        print(sampleSheet)
        return()
    }
    out.dir <- paste0("DBA.", group1, "-", group2)
    if(!file.exists(out.dir)) dir.create(out.dir)
    o.dir <- setwd(out.dir); on.exit(setwd(o.dir))
    dba <- dba(sampleSheet=sampleSheet)
    dba <- dba.count(dba, summits=250)
    dba <- dba.contrast(dba, categories=DBA_CONDITION, minMembers=2)
    dba <- dba.analyze(dba)
    dba.report <- dba.report(dba)
    write.xls(dba.report, file=paste0("DBA.", out.dir, ".report.xls"))
    pdf(paste0(out.dir, ".pdf")); #on.exit(dev.off())
    dba.plotHeatmap(dba)
    #dba.plotVenn(dba)
    dba.plotMA(dba)
    dba.plotVolcano(dba)
    save(dba, dba.report, file="DBA.RData")
} #}}}

#atac.dir <- data.frame(
                       #group=c("BJFF", "KMT2CHET", "KMT2CKO"),
                       #top.path=c(atac.dir0, atac.dir0, atac.dir0),
                       #dir=c("BJFF_atac", "KMT2Chet_atac", "KMT2C-KO_atac")
                       #)
#run.DiffBind("BJFF", "KMT2CKO", atac.dir)
