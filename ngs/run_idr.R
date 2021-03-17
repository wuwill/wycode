source("~/wycode/ngs/find.file.R")
run_idr <- function(peaks, ta, prefix="idr_output", genome="mm10", run=TRUE){ #{{{
    library(data.table)
    tsv <- find.file(file.path("/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data", genome), "*tsv")
    info <- as.data.frame(fread(tsv, header = FALSE))
    chrsz <- info[match("chrsz", info[[1]]), 2]
    blacklist <- info[match("blacklist", info[[1]]), 2]
    py <- "python /opt/apps/labs/gtac/software/atac/miniconda3/envs/encode-chip-seq-pipeline/bin/encode_idr.py"

    cmd <- paste(py, paste0(peaks, collapse = " "), "\\
                  --prefix", prefix, "\\
                  --idr-thresh 0.1 \\
                  --peak-type narrowPeak\\
                  --idr-rank p.value\\
                  --chrsz", chrsz, "\\
                  --blacklist", blacklist, "\\
                   \\
                   --ta", ta)
    sh.file <- paste0(prefix, "_idr.sh")
    cat("#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=ar_minus_input
#SBATCH --mem=24G
#SBATCH --cpus-per-task=3

source /scratch/gtac/software/atac/set.env3.sh
source activate encode-chip-seq-pipeline
cd ", dirname(prefix), "

", cmd,
        file = sh.file, sep="")
    if(run) {
        o.dir <- setwd(dirname(prefix)); on.exit(setwd(o.dir))
        system(paste("sbatch", basename(sh.file)))
    }
    return(sh.file)
} #}}}
run_reproducibility <- function(top.dir, genome="mm10", prefix="idr_output", run = TRUE){ #{{{
    library(data.table)
    tsv <- find.file(file.path("/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data", genome), "*tsv")
    info <- as.data.frame(fread(tsv))
    chrsz <- info[match("chrsz", info[[1]]), 2]
    blacklist <- info[match("blacklist", info[[1]]), 2]
    py <- "python /opt/apps/labs/gtac/software/atac/miniconda3/envs/encode-chip-seq-pipeline/bin/encode_reproducibility_qc.py"
    peaks <- find.file(top.dir, "rep*_rep*.idr*bfilt.narrowPeak.gz")
    peaks.pr <- find.file(top.dir, "rep*pr.idr*bfilt.narrowPeak.gz")
    peaks.ppr <- find.file(top.dir, "ppr.idr*bfilt.narrowPeak.gz")
                 cmd <- paste(py, paste(peaks, collapse = " "), "\\
                              --peaks-pr", paste(peaks.pr, collapse = " "), "\\
                              --peak-ppr", peaks.ppr, "\\
                              --prefix idr \\
                              --peak-type narrowPeak \\
                              \\
                              --chrsz", chrsz)
    exec.dir <- gsub("/IDR/.*$", "/IDR_reproducibility", peaks[1])
    if(!file.exists(exec.dir)) dir.create(exec.dir)
    prefix <- file.path(exec.dir, prefix)

sh.file <- paste0(prefix, "_idr.sh")
cat("#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=ar_minus_input
#SBATCH --mem=24G
#SBATCH --cpus-per-task=3

source /scratch/gtac/software/atac/set.env3.sh
source activate encode-chip-seq-pipeline
cd ", dirname(prefix), "

", cmd,
    file = sh.file, sep="")
    if(run) {
        o.dir <- setwd(dirname(prefix)); on.exit(setwd(o.dir))
        system(paste("sbatch", basename(sh.file)))
    }
    return(sh.file)
} #}}}
run_idr_for_chip <- function(top.dir, genome, run){ #{{{
    library(stringr)
    library(readr)
    peaks <- find.file(top.dir, "*K.narrowPeak.gz", "call-macs2/shard")
    peaks.pr1 <- find.file(top.dir, "*K.narrowPeak.gz", "call-macs2_pr1")
    peaks.pr2 <- find.file(top.dir, "*K.narrowPeak.gz", "call-macs2_pr2")
    peaks.ppr <- find.file(top.dir, "*K.narrowPeak.gz", "call-macs2_ppr")
    peaks.pooled <- find.file(top.dir, "*K.narrowPeak.gz", "call-macs2_pooled")
    ta <- find.file(top.dir, "*tagAlign.gz", "call-bam2ta/shard")
    ta.pooled <- find.file(top.dir, "*tagAlign.gz", "call-pool_ta/exe")
    samples <- str_extract(peaks, "shard-[0-9]+")
    n <- length(peaks)

    exec.dir <- gsub("call-pool_ta/.*$", "", ta.pooled)
    idr.dir <- file.path(exec.dir, "IDR")
    if(!file.exists(idr.dir)) dir.create(idr.dir)

    # pairs of samples
    for(i in 1:(n - 1)) { #{{{
        for(j in (i + 1):n){ #{{{
            run_idr(c(peaks[c(i, j)], peaks.pooled), ta.pooled, paste0(idr.dir, "/rep", i, "_rep", j), genome, run)
        } #}}}
    } #}}}

    # pairs of pr
    for(i in 1:n) 
        run_idr(c(grep(samples[i], peaks.pr1, value=TRUE),
                  grep(samples[i], peaks.pr2, value=TRUE),
                  peaks.pooled),
                ta[i],
                paste0(idr.dir, "/rep", i, "_pr"), genome, run)

    # pprs
        run_idr(c(peaks.ppr, peaks.pooled), ta.pooled, paste0(idr.dir, "/ppr"), genome, run)
} #}}}
run_idr_for_chip("/scratch/gtac/analysis/chip_seq/2898_4_mahajan.072019/merged/ar_merged_minus_input/", "hg19", run=TRUE)
debug(run_reproducibility)
run_reproducibility("/scratch/gtac/analysis/chip_seq/2898_4_mahajan.072019/merged/ar_merged_minus_input/", "hg19")

	#--prefix rep1-pr \
	#s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pr1.pval0.01.300K.narrowPeak.gz
    #s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pr2.pval0.01.300K.narrowPeak.gz
    #s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pval0.01.300K.narrowPeak.gz
	#s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.tagAlign.gz

	#--prefix ppr \
	#s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pr1.pooled.pval0.01.300K.narrowPeak.gz
    #s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pr2.pooled.pval0.01.300K.narrowPeak.gz
    #s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pooled.pval0.01.300K.narrowPeak.gz
	#s_1_1_withindex_sequence.txt_GTAGAGGA.trim.merged.nodup.pooled.tagAlign.gz



