createChipseqInputJson.obsolete <- function(fq1, fq2=NULL, fq1.input=NULL, fq2.input=NULL, genome='mm10',
                                            title="", description="", type="tf",
                                            use_pooled_ctl=FALSE,
                                            file="",
                                            ...){ #{{{
    # by will@wustl.edu at 02/27/2019
    ## other args in ... -
    ## write to json file if file is specified

    args <- list(...)
    library(rjson)
    input1 <- input2 <- input3 <- list()
    if(genome == "mm10") input1$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm10/mm10.tsv"
    if(genome == "hg19") input1$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/hg19/hg19.tsv"
    if(genome == "mm9") input1$chip.genome_tsv <- "/scratch/ref/gtac/reference_sequences/chipseq_pipeline_genome_data/mm9/mm9.tsv"
    input1$chip.paired_end <- is.null(fq2)
    if(is.list(fq1)){ #{{{
        nrep <- length(fq1)
        for(i in 1:nrep){ #{{{
            input2[[paste0("chip.fastqs_rep", i, "_R1")]] <- fq1[[i]]
            if(input2$chip.paired_end)
                input2[[paste0("chip.fastqs_rep", i, "_R2")]] <- fq2[[i]]
        } #}}}
    } else{
        input2[["chip.fastqs_rep1_R1"]] <- fq1
        if(input2$chip.paired_end)
            input2[["chip.fastqs_rep1_R2"]] <- fq2
    } #}}}
    if(is.list(fq1.input)){ #{{{
        nrep <- length(fq1.input)
        for(i in 1:nrep){ #{{{
            input2[[paste0("chip.ctl_fastqs_rep", i, "_R1")]] <- fq1.input[[i]]
            if(input2$chip.paired_end)
                input2[[paste0("chip.ctl_fastqs_rep", i, "_R2")]] <- fq2.input[[i]]
        } #}}}
    } else if(length(fq1.input)>0) {
        input2[["chip.ctl_fastqs_rep1_R1"]] <- fq1.input
        if(input2$chip.paired_end)
            input2[["chip.ctl_fastqs_rep1_R2"]] <- fq2.input
    } #}}}
    input3$chip.title <- title
    input3$chip.description <- description
    input3$chip.pipeline_type <- type
    input3$chip.dup_marker <- "picard"
    input3$chip.mapq_thresh <- 30
    input3$chip.no_dup_removal <- FALSE
    input3$chip.mito_chr_name <- "chrM"
    input3$chip.subsample_reads <- 0
    input3$chip.ctl_subsample_reads <- 0
    input3$chip.xcor_subsample_reads <- 15000000
    input3$chip.keep_irregular_chr_in_bfilt_peak <- FALSE
    input3$chip.always_use_pooled_ctl <- use_pooled_ctl
    input3$chip.ctl_depth_ratio <- 1.2
    input3$chip.macs2_cap_num_peak <- 500000
    input3$chip.pval_thresh <- 0.01
    input3$chip.idr_thresh <- 0.05
    input3$chip.spp_cap_num_peak <- 300000
    input3$chip.bwa_cpu <- 4
    input3$chip.bwa_mem_mb <- 20000
    input3$chip.bwa_time_hr <- 48
    input3$chip.filter_cpu <- 2
    input3$chip.filter_mem_mb <- 20000
    input3$chip.filter_time_hr <- 24
    input3$chip.bam2ta_cpu <- 2
    input3$chip.bam2ta_mem_mb <- 10000
    input3$chip.bam2ta_time_hr <- 6
    input3$chip.spr_mem_mb <- 16000
    input3$chip.fingerprint_cpu <- 2
    input3$chip.fingerprint_mem_mb <- 12000
    input3$chip.fingerprint_time_hr <- 6
    input3$chip.xcor_cpu <- 2
    input3$chip.xcor_mem_mb <- 16000
    input3$chip.xcor_time_hr <- 24
    input3$chip.macs2_mem_mb <- 16000
    input3$chip.macs2_time_hr <- 24
    input3$chip.spp_cpu <- 2
    input3$chip.spp_mem_mb <- 16000
    input3$chip.spp_time_hr <- 72
    input3<- c(input, args)
    ret <- paste0(toJSON(input1, pretty=TRUE, autounbox=TRUE), 
                  toJSON(input2, pretty=TRUE, autounbox=FALSE),
                  toJSON(input3, pretty=TRUE, autounbox=TRUE)
                  )
    ret <- gsub("\n}\n{", ",\n", ret, fixed=TRUE)
    if(!is.null(file) && file!="") cat(ret, file=file)
    return(invisible(ret))
    # use toJSON(input, beatify=TRUE) to convert 
} #}}}

createChipseqInputJson <- function(fq1, fq2=NULL, fq1.input=NULL, fq2.input=NULL, genome='mm10',
                                   title="", description="", type="tf",
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
            input[[paste0("chip.fastqs_rep", i, "_R1")]] <- fq1[[i]]
            if(input$chip.paired_end)
                input[[paste0("chip.fastqs_rep", i, "_R2")]] <- fq2[[i]]
        } #}}}
    } else{
        input[["chip.fastqs_rep1_R1"]] <- fq1
        if(input$chip.paired_end)
            input[["chip.fastqs_rep1_R2"]] <- fq2
    } #}}}
    if(is.list(fq1.input)){ #{{{
        nrep <- length(fq1.input)
        for(i in 1:nrep){ #{{{
            input[[paste0("chip.ctl_fastqs_rep", i, "_R1")]] <- fq1.input[[i]]
            if(input$chip.paired_end)
                input[[paste0("chip.ctl_fastqs_rep", i, "_R2")]] <- fq2.input[[i]]
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
        #SBATCH --mem=", 20*nrep, "G
        #SBATCH --cpus-per-task=", nrep+1, "

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


    sh.file <- createChipseqSlurmSh(json.file, nrep=length(i.sample))
    if(submit) system(paste("sbatch", sh.file))
    return(sh.file)
} #}}}

compileMultiRunResults <- function(top.dir){ #{{{
    result.path <- file.path(top.dir, "cromwell-executions/chip")
    old.dir <- setwd(result.path); on.exit(setwd(old.dir))
    runs <- paste0(system("ls -t", TRUE), "/")
    n <- length(runs)
    if(n<2) return()
    for(i in 1:(n-1)) {
        system(paste0("rsync -a ", runs[i], " ", runs[i+1]))
        system(paste0("rm -rf ", runs[i]))
    }
} #}}}

cpToHTCFdownload <- function(top.dir, dest.dir, files2cp=c('qc.html', "qc.json", "*pooled.fc.signal.bigwig", "*optimal_peak.narrowPeak.*gz", 'overlap.reproducibility.qc')){ #  {{{
                             old.dir <- setwd(top.dir); on.exit(setwd(old.dir))
                             dest.dir <- file.path(dest.dir, basename(top.dir))
                             if(!file.exists(dest.dir)) dir.create(dest.dir)
                             find.files <- function(pattern){ #{{{
                                 system(paste0("find -name '", pattern, "' | grep -v inputs"), TRUE)
                             } #}}}
                             file.copy(unlist(lapply(files2cp, find.files)), dest.dir)
}#}}}
