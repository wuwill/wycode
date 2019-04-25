#!/opt/apps/labs/gtac/opt/R-3.4.1/bin/Rscript
suppressPackageStartupMessages(library(optparse))

# options --------------------------------------------------------------------------------#
option_list <- list(
make_option(c("-r", "--run"), action="store", help = "run_lanes in the form of 1234_12; for multiple runs, quote with ''."),
make_option(c("-g", "--genome"), action = "store", help = "Reference genome: hg19, mm9, or mm10."),
make_option(c("-c", "--customer"), action = "store", help = "Customer."),
make_option(c("-k", "--samplekey"), action = "store", default=NULL, help = "Sample key file. Default: look within GTAC sample key directory"),
make_option(c("-G", "--group"), action = "store", default=NULL, help = "Group information file. Default: 'samples.txt' within the project directory.
                    Tab delimited file with either 2 or 4 columns.
                    - 2 columns: sample name, sample group;
                    - 4 columns: sample key with additional 4th column for sample group.
                    Alternatively, use tab delimited file with header - name, group, fq1 (and/or) fq2. This allows using custom fastq files as input."),
make_option(c("-C", "--contrast"), action = "store", default=NULL, help = "Contrast groups file. Default: 'contrast.txt' within the project directory.
                    Tab delimited file with 2 columns: sample group, control group."),
make_option(c("--notsubmit"), action = "store_true", default=FALSE, help = "When specified, creates script files without submitting jobs to SLURM."),
make_option(c("-d", "--directory"), action = "store", default=NULL, help = "Work directory within '/scratch/gtac/analysis/chipseq'. Default: [run_lane_customer]."),
make_option(c("-p", "--postresults"), action = "store_true", default=FALSE, help = "Post analysis results rather than perform analysis."),
make_option(c("-o", "--overwrite"), action = "store_true", default=FALSE, help = "- For posting resutls: when specified, overwrite posted result files."),
make_option(c("-u", "--url"), action = "store", default=NULL, help = "- For posting resutls: GTAC_RUN_DATA URL. Defult: generate URL using HTCF 'serve'. Make sure to 'module load HTCF'.")
)
opt_parser <- OptionParser(usage = "%prog [options]\n\nExamples:\n%prog -r 1234_12 -c customer -g hg19 --notsubmit\n%prog -r 1234_12 -c customer -g hg19 --postresults\n\n>>> For more complex situations, make sure to specify sample groups and contrasting groups by creating samples.txt and contrast.txt within the project directory before running this command (see options below) !!!", option_list=option_list)
opt <- parse_args(opt_parser)
if(is.null(opt$run) && is.null(opt$directory)) {
    print_help(opt_parser)
    stop("More options required for running the script.")
}
submit <- !opt$notsubmit

runs <- strsplit(opt$run, " +")[[1]]
genome <- opt$genome
customer <- tolower(opt$customer)
POST_RESULTS <- opt$postresults # MAKE SURE to 'module load HTCF' before running R script to post results, so that URL will be created for the posted data.
sample.key.file <- if(is.null(opt$samplekey)) opt$samplekey else normalizePath(opt$samplekey)
sample.group.file <- if(is.null(opt$group)) opt$group else normalizePath(opt$group)
comparison.file <- if(is.null(opt$contrast)) opt$contrast else normalizePath(opt$contrast)
work.dir <- opt$directory
cat("Specified options:\n")
print(opt)

# deparse multiple run lanes --------------------------------------------------------------------------------#
if(length(runs)>1){
    run_lanes <- strsplit(runs, "_")
    runs <- sapply(run_lanes, function(x) x[1])
    lanes <- sapply(run_lanes, function(x) x[2])
    if(any(nchar(lanes)>1)) {
        lanes <- strsplit(lanes, "")
        runs <- rep(runs, sapply(lanes, length))
        lanes <- unlist(lanes)
    }
    run_lanes <- cbind(runs, lanes)
    run_lanes <- run_lanes[order(runs, lanes),,drop=FALSE]
    RUNS <- by(run_lanes, run_lanes[,1], function(x) paste(x[1, 1], paste(x[,2], collapse=""), sep="_"))
    RUNS <- paste(RUNS, collapse="_")
    runs <- apply(run_lanes, 1, paste, collapse="_")
} else {
    RUNS <- runs
}


work.dir <- if(is.null(work.dir)) file.path("/scratch/gtac/analysis/chip_seq", paste0(RUNS, "_", customer)) else
    file.path("/scratch/gtac/analysis/chip_seq", work.dir)
cat("\nWork directory:", work.dir, "\n")
if(!file.exists(work.dir)) {
    cat(" - Work directory does not exist. Creating work directory.\n")
    dir.create(work.dir)
}
setwd(work.dir)
demux.dir <- "/scratch/gtac/analysis/demux"
source("/scratch/gtac/software/atac/chipseq2.launch.func.r")
library(readr)

#--------------------------------------------------------------------------------#
# CREATE ANALYSIS PROJECT and RUN the PIPELINE
if(!POST_RESULTS) {  #{{{
    # look for sample.group.file & comparison.file
    if(is.null(sample.group.file) && file.exists("samples.txt")) sample.group.file <- "samples.txt"
    if(is.null(comparison.file) && file.exists("contrasts.txt")) comparison.file <- "contrasts.txt"

    if(!is.null(sample.group.file)){
        cat("Reading sample group info from", sample.group.file, "\n")
        group.info <- read_tsv(sample.group.file, col_names=FALSE)
        if(all(c("name", "group") %in% unlist(group.info[1,]))) {
            cat(" - The sample group info file has HEADERS. Use it to annotate sample and their groups.\n")
            info <- read_tsv(sample.group.file, col_names=TRUE)
    } else
            {
                if(ncol(group.info)==2) {
                    cat(" -The sample group info file has 2 columns.\n")
                    if(is.null(sample.key.file)) sample.key.file <- Sys.glob(paste0("/scratch/ref/gtac/samplekeys/", runs[1], "*_", customer, "*"))
                    cat("Read sample info from", sample.key.file, "\n")
                    info <- read_tsv(sample.key.file, col_names=FALSE)
                    colnames(info) <- if(ncol(info)==4)
                        c("id", "index", "name", "group") else
                            c("id", "index", "name")
                    info$group <- group.info$X2[match(info$name, group.info$X1)]
                } else if(ncol(group.info)==4) {
                    cat(" - The sample group info file has 4 columns. Use this table to annotate sample and their groups.\n")
                    colnames(group.info) <- c("id", "index", "name", "group")
                    info <- group.info
                } else {
                    cat(" - Use the sample group file as sample key and treat each sample as group by itself.\n")
                    group.info <- group.info[,1:3,drop=FALSE]
                    colnames(group.info) <- c("id", "index", "name")
                    group.info$group <- group.info$name
                    info <- group.info
                }
            }
    } else{
        cat("No sample group info file has been specified/found. Use information from sample key file.\n")
        if(is.null(sample.key.file)) sample.key.file <- Sys.glob(paste0("/scratch/ref/gtac/samplekeys/", runs[1], "*_", customer, "*"))
        cat(" - Read sample info from", sample.key.file, "\n")
        info <- read_tsv(sample.key.file, col_names=FALSE)
        colnames(info) <- if(ncol(info)==4)
            c("id", "index", "name", "group") else
                c("id", "index", "name")

    }
    if(is.null(info$group)) {
        cat("No group information provided! Treat each sample as a group by itself!\n")
        info$group <- info$name
    }

    if(!is.null(comparison.file)){
        cat("Read contrasting groups from", comparison.file, "\n")
        contrastingGroups <- read_tsv(comparison.file, col_names=FALSE)
    } else {
        cat("No contrasting information provided! Try figuring out contrasting groups by comparing to 'input' and 'igg'!\n")
        info$group <- tolower(info$group)
        control <- info$group[info$group %in% c("input", "igg")]
        contrastingGroups <- cbind(rep(setdiff(info$group, control), each=length(control)), control)
        if(nrow(contrastingGroups)<1) stop(" - Failed to figure out contrasting groups. Fix this by providing a contrast groups file using '-C'.")
    }

    cat("Sample information:\n")
    print(as.data.frame(info))
    cat("\nContrast groups:\n")
    print(contrastingGroups)

    for(i in 1:nrow(contrastingGroups)){ #{{{
        task.name <- paste0(contrastingGroups[i,1], "_minus_", contrastingGroups[i,2])
        createChipseqFromFq(contrastingGroups[i,1], contrastingGroups[i,2],
                            info, genome, file.path(demux.dir, runs), type="histone",
                            file=task.name, submit=submit)
    } #}}}
} #}}}

#--------------------------------------------------------------------------------#
# post anlysis results
if(POST_RESULTS){
    cat("Make sure 'module load htcf' before posting results!\n")
    library(R.utils)
    overwrite <- opt$overwrite
    top.url <- opt$url
    post.dir <- file.path("/scratch/gtac/GTAC_RUN_DATA/", paste0(capitalize(customer), '_', RUNS))
    postResults(dest.dir=post.dir, overwrite=overwrite)
    postDemux(runs, dest.dir=post.dir, overwrite=overwrite)
    if(is.null(top.url)) top.url <- serveData(post.dir)
    cat("Data available at:\n", top.url, "\n\n")
    cat("QC reports available at:", "\n")
    cat(getQcUrl(post.dir, top.url=top.url), sep="\n")
    cat("\n")
    cat("Siganal tracks could be viewed on the WU epigenome browser from", "\n")
    getWuBrowserLink(post.dir, top.url=top.url, genome=genome)
    cat("\n")
}
