createChipseqInputJson <- function(fq1, fq2=NULL, fq1.input=NULL, fq2.input=NULL, genome='mm10',
                                   title=basename(file), description=title, type="tf",
                                   use_pooled_ctl=FALSE,
                                   file="",
                                   ...){ #{{{
    # by will@wustl.edu at 02/27/2019
    ## other args in ... -
    ## write to json file if file is specified

    args <- list(...)
    library(rjson)
    input <- list()
    input1 <- input2 <- input3 <- list()
    if(genome == "mm10") input$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm10/mm10.tsv"
    if(genome == "hg19") input$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/hg19/hg19.tsv"
    if(genome == "mm9") input$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm9/mm9.tsv"
    input$chip.paired_end <- !is.null(fq2)
    if(is.list(fq1)){ #{{{
        nrep <- length(fq1)
        for(i in 1:nrep){ #{{{
            input[[paste0("chip.fastqs_rep", i, "_R1")]] <- if(nrep>1) fq1[[i]] else fq1
            if(input$chip.paired_end)
                input[[paste0("chip.fastqs_rep", i, "_R2")]] <- if(nrep>1) fq2[[i]] else fq2
        } #}}}
    } else{
        input[["chip.fastqs_rep1_R1"]] <- fq1
        if(input$chip.paired_end)
            input[["chip.fastqs_rep1_R2"]] <- fq2
    } #}}}
    if(is.list(fq1.input)){ #{{{
        nrep <- length(fq1.input)
        for(i in 1:nrep){ #{{{
            input[[paste0("chip.ctl_fastqs_rep", i, "_R1")]] <- if(nrep>1) fq1.input[[i]] else fq1.input
            if(input$chip.paired_end)
                input[[paste0("chip.ctl_fastqs_rep", i, "_R2")]] <- if(nrep>1) fq2.input[[i]] else fq2.input
        } #}}}
    } else if(length(fq1.input)>0) {
        input[["chip.ctl_fastqs_rep1_R1"]] <- fq1.input
        if(input$chip.paired_end)
            input[["chip.ctl_fastqs_rep1_R2"]] <- fq2.input
    } #}}}
    input$chip.title <- title
    input$chip.description <- description
    input$chip.pipeline_type <- type
    #input$chip.dup_marker <- "picard"
    input$chip.always_use_pooled_ctl <- use_pooled_ctl
    input<- c(input, args)
    ret <- toJSON(input, indent=4)
    if(!is.null(file) && file!="") cat(ret, file=file)
    return(invisible(ret))
    # use toJSON(input, beatify=TRUE) to convert 
} #}}}

createChipseqSlurmSh <- function(input.json, nrep=2){ #{{{
    job.name <- gsub(".json$", "", basename(input.json))
    sh.file <- gsub(".json$", ".sh", input.json)
    parent.dir <- dirname(normalizePath(input.json))

    cat("#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=", job.name, "
#SBATCH --mem=", min(max(12*nrep, 24), 64), "G
#SBATCH --cpus-per-task=", min(nrep+1, 5), "

cd ", parent.dir, "
module load java || true
source /scratch/gtac/software/atac/set.env.sh
source activate encode-chip-seq-pipeline

java -jar -Dconfig.file=/scratch/gtac/software/atac/chip-seq-pipeline2/backends/backend.conf \\
-Dbackend.providers.Local.config.concurrent-job-limit=", nrep, " \\
/scratch/gtac/software/atac/cromwell-34.jar run /scratch/gtac/software/atac/chip-seq-pipeline2/chip.wyang.wdl -i ", basename(input.json), " -m metadata.json
", file=sh.file, sep="")
    return(sh.file)
} #}}}

findFq <- function(index, dirs, r1.pattern="_1_withindex"){ #{{{
    fqs <- lapply(dirs,function(x) Sys.glob(file.path(x, paste0("*", index, "*q.gz"))))
    if(length(fqs[[1]])==1){ #{{{
        return(list(fq1=unlist(fqs)))
    } #}}}
    fq1 <- lapply(fqs, grep, pattern=r1.pattern, value=TRUE)
    fq2 <- mapply(setdiff, fqs, fq1)
    return(list(fq1=unlist(fq1), fq2=unlist(fq2)))
} #}}}

createChipseqFromFq <- function(group, input.group, sample.info, genome, fq.dirs, type="tf", suffix="", file=paste(group, input.group, sep="_x_", ifelse(suffix=="", "", paste0("_", suffix))), submit=TRUE, ...){ #{{{
    ## sample.info columns
    #   - manditory: group, name
    #   - optional:  (fq1 | index)
    i.sample <- sample.info$group %in% group
    i.input <- sample.info$group %in% input.group
    if(is.null(sample.info$fq1)){  #{{{ ## find fq files based on indexes
        sample.indexes <- sample.info$index[i.sample]
        input.indexes <- sample.info$index[i.input]
        fqs <- lapply(sample.indexes, findFq, dirs=fq.dirs)
        fqs.input <- lapply(input.indexes, findFq, dirs=fq.dirs)
        fq1 <- lapply(fqs, function(x) x$fq1)
        fq2 <- lapply(fqs, function(x) x$fq2)
        fq1.input <- lapply(fqs.input, function(x) x$fq1)
        fq2.input <- lapply(fqs.input, function(x) x$fq2)
        if(length(unlist(fq2))<1)
            fq2 <- fq2.input <- NULL
        #}}}
    } else { #  {{{ ## get fq files from columns fq1, fq2
        get.parsed.fq <- function(i, forward=TRUE){ #{{{
            fqs <- if(forward) sample.info$fq1[i] else sample.info$fq2[i]
            if(any(grepl(";", fqs))) return(lapply(strsplit(fqs, ";"), gsub, pattern=" +", replacement=""))
            if(any(grepl(",", fqs))) return(lapply(strsplit(fqs, ","), gsub, pattern=" +", replacement=""))
            if(any(grepl("\t", fqs))) return(lapply(strsplit(fqs, "\t"), gsub, pattern=" +", replacement=""))
            if(any(grepl(" ", fqs))) return(lapply(strsplit(fqs, " +"), function(x) x))
            return(lapply(strsplit(fqs, " +"), function(x) x))
        } #}}}
        fq1 <- lapply(which(i.sample), get.parsed.fq)
        fq1.input <- lapply(which(i.input), get.parsed.fq)
        if(is.null(sample.info$fq2)){
            fq2 <- fq2.input <- NULL
        } else {
            fq2 <- lapply(which(i.sample), get.parsed.fq, forward=FALSE)
            fq2.input <- lapply(which(i.input), get.parsed.fq, forward=FALSE)
        }
    } #}}}

    dir.create(file)
    json.file <- paste0(file, "/", file, ".json")

    createChipseqInputJson(fq1, fq2, fq1.input, fq2.input, genome=genome,
                           type=type,
                           use_pooled_ctl=TRUE,
                           file=json.file, ...)


    sh.file <- createChipseqSlurmSh(json.file, nrep=sum(i.sample))
    if(submit) system(paste("sbatch", sh.file))
    return(sh.file)
} #}}}

compileMultiResumedRuns <- function(top.dir){ #{{{
    result.path <- Sys.glob(file.path(top.dir, "cromwell-executions/*"))
    old.dir <- setwd(result.path); on.exit(setwd(old.dir))
    runs <- paste0(system("ls -t", TRUE), "/")
    n <- length(runs)
    if(n<2) return()
    for(i in 1:(n-1)) {
        system(paste0("rsync -a ", runs[i], " ", runs[i+1]))
        system(paste0("rm -rf ", runs[i]))
    }
} #}}}

cpToHTCFdownload <- function(top.dir, dest.dir, files2cp=c('qc.html', "qc.json", "*pooled.fc.signal.bigwig", "*optimal*.narrowPeak.*gz", 'overlap.reproducibility.qc')){ #  {{{
                             old.dir <- setwd(top.dir); on.exit(setwd(old.dir))
                             dest.dir <- file.path(dest.dir, basename(top.dir))
                             if(!file.exists(dest.dir)) dir.create(dest.dir)
                             find.files <- function(pattern){ #{{{
                                 system(paste0("find -name '", pattern, "' | grep -v inputs"), TRUE)
                             } #}}}
                             file.copy(unlist(lapply(files2cp, find.files)), dest.dir)
}#}}}
postResults <- function(top.dir=NULL, dest.dir, overwrite=FALSE, atac=FALSE){ #{{{
    if(is.null(top.dir)){ #{{{
        top.dirs <- list.files(pattern=ifelse(atac, "_atac$", "_minus_"))
        print(top.dirs)
        sapply(top.dirs, postResults, dest.dir=dest.dir, overwrite=overwrite)
        return()
    } #}}}
    if(!file.exists(dest.dir)) dir.create(dest.dir)
    result.dir <- file.path(dest.dir, basename(top.dir))
    if(!overwrite && file.exists(result.dir)) {
        cat("Destination directory", result.dir, "exists. \n\tSkip posting data. Set the overwrite option to over write results!")
        return()
    }
    compileMultiResumedRuns(top.dir)
    old.dir <- setwd(top.dir); on.exit(setwd(old.dir))
    if(!file.exists(result.dir)) dir.create(result.dir)
    system(paste0("cp -r cromwell-executions/*/*/call-* ", result.dir))
    system(paste0("cp *.json ", result.dir))
    return(invisible(TRUE))
} #}}}
postDemux <- function(runs, dest.dir, overwrite=FALSE){ #{{{
    if(!file.exists(dest.dir)) dir.create(dest.dir)
    dest.dir <- file.path(dest.dir, "demux")
    if(!overwrite && file.exists(dest.dir)) {
        cat("Destination directory", dest.dir, "exists. \n\tSkip posting data. Set the overwrite option to over write results!")
        return()
    }
    if(!file.exists(dest.dir)) dir.create(dest.dir)
    if(file.exists(file.path("/scratch/gtac/analysis/demux/", runs[1])))
        sapply(runs, function(run) system(paste0("cp -r /scratch/gtac/analysis/demux/", run, " ", dest.dir)))
} #}}}

getQcUrl <- function(dest.dir, top.url=NULL){ #{{{
    qc.html <- Sys.glob(file.path(dest.dir, "*/call-qc_report/execution/qc.html"))
    qc.idr <- Sys.glob(file.path(dest.dir, "*/call-reproducibility_overlap/execution/overlap.reproducibility.qc"))
    ret <- c(qc.html, qc.idr)
    cat(gsub(dest.dir, top.url, ret, fix=TRUE), sep="\n")
    return(invisible(ret))
} #}}}

find.file <- function(path, pattern, pattern2=NULL, vpattern=NULL){ #{{{
    cmd <- paste0("find -H ", path, " -name '", pattern, "'")
    if(!is.null(pattern2)) cmd <- paste0(cmd, " | grep '", pattern2, "'")
    if(!is.null(vpattern)) cmd <- paste0(cmd, " | grep -v '", vpattern, "'")
    return(system(cmd, TRUE))
} #}}}
serveData <- function(dest.dir){ #{{{
    system("module load htcf") # module has to be loaded before starting R session to serve directories
    o.dir <- setwd("/scratch/gtac/GTAC_RUN_DATA"); on.exit(setwd(o.dir))
    system(paste0("chmod -R +r ", dest.dir))
    link <- system(paste("serve", basename(dest.dir)), TRUE)
    out <- c("/scratch/gtac/GTAC_RUN_DATA", basename(dest.dir), "", "", link, format(Sys.Date(), "%m/%d/%y"))
    out.file <- "SERVER_URLs"
    cat(out, file=out.file, sep="\t", append=TRUE)
    cat("\n", file=out.file, append=TRUE)
    print(link)
    return(invisible(link))
} #}}}
getWuBrowserLink <- function(dest.dir, top.url=NULL, genome="hg19"){ #{{{
    dest.dir <- gsub("/+$", "", dest.dir)
    bigwig0 <- find.file(dest.dir, '*.fc.signal.bigwig', "call-macs2/shard-0")
    chip.subfolders <- gsub("/cromwell-ex.*$", "", bigwig0)
    bigwig <- sapply(seq_along(chip.subfolders), function(i){
                     ret <- suppressWarnings(ret <- find.file(chip.subfolders[i], "*.fc.signal.*", "call-macs2.*pooled/"))
                     if(length(ret)<1) ret <- bigwig0[i]
                     return(ret)})
    #bigwig <- find.file(dest.dir, '*.fc.signal.bigwig', "call-macs2.*pooled/")
    json.file <- file.path(dest.dir, "pool_fc.json")
    library(rjson)
    get.track.json <- function(bigwig){ #{{{
        name <- strsplit(gsub(paste0(dest.dir, "/"), "", bigwig), "/")[[1]][1]
        toJSON(list(type="bigwig",
             name=name,
             url=gsub(dest.dir, top.url, bigwig, fix=TRUE),
             mode=1,
             qtc=list(height=30, summeth=2, smooth=3, pr=255, pg=20, pb=147, thtype=1, thmin=2, thmax=40)
             ))
    } #}}}
    json <- sapply(bigwig, get.track.json)
    json <- c(json,
              toJSON(list(type="native_track",
                   list=list(list(name="refGene", mode="full"))
                   ))
              )
    cat("[", paste(json, collapse=",\n"), "]", file=json.file, sep="\n")
    url <- paste0("http://epigenomegateway.wustl.edu/browser/?genome=", genome, "&datahub=", gsub(dest.dir, top.url, json.file, fix=TRUE))
    cat(url)
    return(url)
} #}}}
getWuBrowserLink4bw <- function(bigwig, track.names=gsub(".fc.signal.b.*$", "", basename(bigwig)), json.file='', dest.dir, top.url=NULL, genome="hg19"){ #{{{
    bigwig <- normalizePath(bigwig)
    json.file <- if(is.null(json.file) || json.file %in% c("", NA)) file.path(dest.dir, "pool_fc.json") else file.path(dest.dir, json.file)
    library(rjson)
    # check if bigwig files are within dest.dir; if not create a sub folder called 'signals', and copy bigwig files
    if(any(!grepl(dest.dir, bigwig, fix=TRUE))) { #{{{
        bigwig.dir <- file.path(dest.dir, "signals")
        if(!file.exists(bigwig.dir)) dir.create(bigwig.dir)
        to.files <- file.path(bigwig.dir, paste0(track.names, '.fc.signal.bigwig'))
        file.copy(bigwig, to.files)
        bigwig <- to.files
    } #}}}
    get.track.json <- function(bigwig, track.names){ #{{{
        toJSON(list(type="bigwig",
             name=track.names,
             url=gsub(dest.dir, top.url, bigwig, fix=TRUE),
             mode=1,
             qtc=list(height=30, summeth=2, smooth=3, pr=255, pg=20, pb=147, thtype=1, thmin=2, thmax=40)
             ))
    } #}}}
    #json <- sapply(bigwig, get.track.json, )
    json <- mapply(get.track.json, bigwig, track.names)
    json <- c(json,
              toJSON(list(type="native_track",
                   list=list(list(name="refGene", mode="full"))
                   ))
              )
    cat("[", paste(json, collapse=",\n"), "]", file=json.file, sep="\n")
    url <- paste0("http://epigenomegateway.wustl.edu/browser/?genome=", genome, "&datahub=", gsub(dest.dir, top.url, json.file, fix=TRUE))
    cat(url, "\n")
    return(url)
} #}}}

### epic2
epic2 <- "/opt/apps/labs/gtac/software/atac/miniconda_/bin/epic2"
tabixBed <- function(bed, out=paste0(bed, ".gz"), run=FALSE) { #{{{
    ret <- paste0("(head -n 1 ", bed, " && tail -n +2 ", bed, " | sort -k1,1 -k2,2n ) | /opt/apps/labs/gtac/software/atac/miniconda_/bin/bgzip > ", out, "\n",
                  "/opt/apps/labs/gtac/software/atac/miniconda_/bin/tabix -p bed ", out, "\n")
    if(run) system(ret)
    return(invisible(ret))
} #}}}
createEpic2SlurmSh <- function(chip.dir, genome="hg19", submit=FALSE){ #{{{
    chip.dir <- normalizePath(chip.dir)
    work.dir <- dirname(find.file(chip.dir, "call-macs2"))
    if(length(work.dir)==0 || work.dir=="") {
        cat("Run chipseq2 pipeline before running epic2 for", chip.dir, "\n")
        return()
    }
    if(length(find.file(chip.dir, "call-pool_ta"))>0) {
        ppr1.ta <- Sys.glob(file.path(work.dir, "call-pool_ta_pr1/execution/*.pooled.tagAlign.gz"))
        ppr2.ta <- Sys.glob(file.path(work.dir, "call-pool_ta_pr2/execution/*.pooled.tagAlign.gz"))
        ctl.ta <- Sys.glob(file.path(work.dir, "call-pool_ta_ctl/execution/*.pooled.tagAlign.gz"))
        pooled.ta <- Sys.glob(file.path(work.dir, "call-pool_ta/execution/*.pooled.tagAlign.gz"))
    } else {
        ppr1.ta <- Sys.glob(file.path(work.dir, "call-*/execution/*.pr1.tagAlign.gz"))
        ppr2.ta <- Sys.glob(file.path(work.dir, "call-*/execution/*.pr2.tagAlign.gz"))
        ctl.ta <- Sys.glob(file.path(work.dir, "call-pool_ta_ctl/execution/*.tagAlign.gz"))
        pooled.ta <- Sys.glob(file.path(work.dir, "call-bam2ta/shard-0/execution/*.tagAlign.gz"))
    }
    sub.jobs <- list("call-macs2_ppr1"=c(ppr1.ta, ctl.ta),
                     "call-macs2_ppr2"=c(ppr2.ta, ctl.ta),
                     "call-macs2_pooled"=c(pooled.ta, ctl.ta))
    mk.epic2.script <- function(sub.job){ #{{{
        ta <- sub.jobs[[sub.job]][1]
        ctl.ta <- setdiff(sub.jobs[[sub.job]], ta)
        result <- file.path(sub.job, "execution/epic2.results.txt")
        ret <- paste0(epic2, " \\\n -t ", ta, " \\\n ",
                      ifelse(length(ctl.ta)==0 || ctl.ta %in% c("", NA), "", paste0("-c ", ctl.ta, " \\\n ")),
                      "--genome ", genome, " \\\n ",
                      "-o ", result, "\n",
                      tabixBed(result),
                      "touch ", file.path(dirname(result), "epic2.done"), "\n"
                      )
        return(ret)
    } #}}}

    sh.file <- file.path(work.dir, "epic2.sh")
    job.name <- paste0("epic_", basename(chip.dir))
    o.dir <- setwd(work.dir); on.exit(setwd(o.dir))

    cat("#!/bin/bash
#SBATCH --job-name=", job.name, "
#SBATCH --mem=15G
source /scratch/gtac/software/atac/set.env.sh
cd ", work.dir, "
rm -f call-macs2_pooled/execution/epic2.done

", sapply(names(sub.jobs), mk.epic2.script), "

cd ", work.dir, "
## wait for epic jobs to finish
#while [! -f call-macs2_pooled/execution/epic2.done]
#do
    #sleep 120
#done 

# get naive overlaps of epic2 peaks
/opt/apps/labs/gtac/opt/R-3.4.1/bin/Rscript /scratch/gtac/software/atac/peak_naive_overlaps.R call-macs2_pooled/execution/epic2.results.txt call-macs2_ppr1/execution/epic2.results.txt call-macs2_ppr2/execution/epic2.results.txt call-macs2_pooled/execution/epic2.naive_overlaps.txt
", tabixBed("call-macs2_pooled/execution/epic2.naive_overlaps.txt"), "
", file=sh.file, sep="")
    if(submit) system(paste0("sbatch ", sh.file))
    return(sh.file)
} #}}}

### ATAC
createAtacInputJson <- function(fq1, fq2=NULL, genome='mm10',
                                   title=basename(file), description=title, type="atac",
                                   file="",
                                   ...){ #{{{
    # by will@wustl.edu at 05/07/2019
    ## other args in ... -
    ## write to json file if file is specified

    args <- list(...)
    library(rjson)
    input <- list()
    input1 <- input2 <- input3 <- list()
    if(genome == "mm10") input$atac.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm10/mm10.tsv"
    if(genome == "hg19") input$atac.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/hg19/hg19.tsv"
    if(genome == "mm9") input$atac.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm9/mm9.tsv"
    input$atac.paired_end <- !is.null(fq2)
    if(is.list(fq1)){ #{{{
        nrep <- length(fq1)
        for(i in 1:nrep){ #{{{
            input[[paste0("atac.fastqs_rep", i, "_R1")]] <- if(nrep>1) fq1[[i]] else fq1
            if(input$atac.paired_end)
                input[[paste0("atac.fastqs_rep", i, "_R2")]] <- if(nrep>1) fq2[[i]] else fq2
        } #}}}
    } else{
        input[["atac.fastqs_rep1_R1"]] <- fq1
        if(input$atac.paired_end)
            input[["atac.fastqs_rep1_R2"]] <- fq2
    } #}}}
    input$atac.auto_detect_adapter <- TRUE
    input$atac.title <- title
    input$atac.description <- description
    input$atac.pipeline_type <- type
    input$"atac.cutadapt_param" <- "-e 0.1 -m 5"
    input$"atac.enable_count_signal_track" <- TRUE
    input$"atac.multimapping" <- 0
    input$"atac.bowtie2_param_pe" <- "-X2000 --mm --local"
    input$"atac.bowtie2_param_se" <- "--local"
    #"atac.dup_marker" : "picard",
    input$"atac.mapq_thresh" <- 30
    input$"atac.no_dup_removal" <- FALSE
    input$"atac.mito_chr_name" <- "chrM"
    input$"atac.regex_filter_reads" <- "chrM"
    input$"atac.subsample_reads" <- 0
    input$"atac.enable_xcor" <- TRUE
    input$"atac.xcor_subsample_reads" <- 25000000
    input$"atac.keep_irregular_chr_in_bfilt_peak" <- FALSE
    input$"atac.cap_num_peak" <- 300000
    input$"atac.pval_thresh" <- 0.01
    input$"atac.smooth_win" <- 150
    input$"atac.enable_idr" <- TRUE
    input$"atac.idr_thresh" <- 0.1
    input$"atac.disable_ataqc" <- FALSE
    input$"atac.bowtie2_cpu" <- 6
    input$"atac.bowtie2_mem_mb" <- 24000
    input$"atac.bowtie2_time_hr" <- 48
    input$"atac.filter_cpu" <- 2
    input$"atac.filter_mem_mb" <- 20000
    input$"atac.filter_time_hr" <- 24
    input$"atac.bam2ta_cpu" <- 2
    input$"atac.bam2ta_mem_mb" <- 10000
    input$"atac.spr_mem_mb" <- 16000
    input$"atac.xcor_cpu" <- 2
    input$"atac.xcor_mem_mb" <- 16000
    input$"atac.xcor_time_hr" <- 6
    input$"atac.macs2_mem_mb" <- 16000
    input$"atac.macs2_time_hr" <- 24
    input$"atac.ataqc_mem_mb" <- 16000
    input$"atac.ataqc_mem_java_mb" <- 15000
    input$"atac.ataqc_time_hr" <- 24
    #input$input$chip.always_use_pooled_ctl <- use_pooled_ctl
    input<- c(input, args)
    ret <- toJSON(input, indent=4)
    if(!is.null(file) && file!="") cat(ret, file=file)
    return(invisible(ret))
    # use toJSON(input, beatify=TRUE) to convert 
} #}}}

createAtacSlurmSh <- function(input.json, nrep=2){ #{{{
    job.name <- gsub(".json$", "", basename(input.json))
    sh.file <- gsub(".json$", ".sh", input.json)
    parent.dir <- dirname(normalizePath(input.json))

    cat("#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=", job.name, "
#SBATCH --mem=", min(max(12*nrep, 24), 64), "G
#SBATCH --cpus-per-task=", min(nrep+1, 5), "

cd ", parent.dir, "
module load java || true
source /scratch/gtac/software/atac/set.env.mc3.sh
source activate encode-atac-seq-pipeline

java -jar -Dconfig.file=/scratch/gtac/software/atac/atac-seq-pipeline/backends/backend.conf \\
-Dbackend.providers.Local.config.concurrent-job-limit=", nrep, " \\
/scratch/gtac/software/atac/cromwell-34.jar run /scratch/gtac/software/atac/atac-seq-pipeline/atac.wyang.wdl -i ", basename(input.json), " -m metadata.json
", file=sh.file, sep="")
    return(sh.file)
} #}}}
createAtacFromFq <- function(group, sample.info, genome, fq.dirs, type="tf", suffix="", file=paste0(group, ifelse(suffix=="", "", paste0("_", suffix))), submit=TRUE, ...){ #{{{
    ## sample.info columns
    #   - manditory: group, name
    #   - optional:  (fq1 | index)
    i.sample <- sample.info$group %in% group
    if(is.null(sample.info$fq1)){  #{{{ ## find fq files based on indexes
        sample.indexes <- sample.info$index[i.sample]
        fqs <- lapply(sample.indexes, findFq, dirs=fq.dirs)
        fq1 <- lapply(fqs, function(x) x$fq1)
        fq2 <- lapply(fqs, function(x) x$fq2)
        if(length(unlist(fq2))<1)
            fq2 <- NULL
        #}}}
    } else { #  {{{ ## get fq files from columns fq1, fq2
        get.parsed.fq <- function(i, forward=TRUE){ #{{{
            fqs <- if(forward) sample.info$fq1[i] else sample.info$fq2[i]
            if(any(grepl(";", fqs))) return(lapply(strsplit(fqs, ";"), gsub, pattern=" +", replacement=""))
            if(any(grepl(",", fqs))) return(lapply(strsplit(fqs, ","), gsub, pattern=" +", replacement=""))
            if(any(grepl("\t", fqs))) return(lapply(strsplit(fqs, "\t"), gsub, pattern=" +", replacement=""))
            if(any(grepl(" ", fqs))) return(lapply(strsplit(fqs, " +"), function(x) x))
            return(lapply(strsplit(fqs, " +"), function(x) x))
        } #}}}
        fq1 <- lapply(which(i.sample), get.parsed.fq)
        if(is.null(sample.info$fq2)){
            fq2 <- NULL
        } else {
            fq2 <- lapply(which(i.sample), get.parsed.fq, forward=FALSE)
        }
    } #}}}

    dir.create(file)
    json.file <- paste0(file, "/", file, ".json")

    createAtacInputJson(fq1, fq2, genome=genome,
                           type=type,
                           file=json.file, ...)


    sh.file <- createAtacSlurmSh(json.file, nrep=sum(i.sample))
    if(submit) system(paste("sbatch", sh.file))
    return(sh.file)
} #}}}
