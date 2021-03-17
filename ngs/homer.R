REF_DIR <- if(file.exists("/scratch/ref")) "/scratch/ref/gtac/reference_sequences" else "~/local"
sys.cmd <- function(...) system(paste0(...))
get.ref.fa <- function(genome=NULL){ #{{{
    if(exists("GENOME") && is.null(genome)) genome <- GENOME
    chip_ref_dir <- file.path(REF_DIR, "chipseq_pipeline_genome_data")
    fa <- if("mm9" %in% genome) file.path(chip_ref_dir, "mm9/mm9.fa") else
        if("mm10" %in% genome) file.path(chip_ref_dir, "mm10/mm10_no_alt_analysis_set_ENCODE.fasta") else
        if("hg19" %in% genome) file.path(chip_ref_dir, "hg19/male.hg19.fa")
    return(fa)
} #}}}
write.peak4homer <- function(peak, file="", rm.dup=FALSE){ #{{{
    seqnm <- if("GRanges" %in% class(peak)) seqnames(peak) else peak$seqnames
    start <- if("GRanges" %in% class(peak)) start(peak) else peak$start
    if(rm.dup){ #{{{
       peak$loc <- paste(seqnm, start, sep="_")
       peak <- peak[!duplicated(peak$loc)]
       seqnm <- if("GRanges" %in% class(peak)) seqnames(peak) else peak$seqnames
       start <- if("GRanges" %in% class(peak)) start(peak) else peak$start
    } #}}}
    peak$name <- paste(peak$symbol, seqnm, start, sep="_")
    if(file=="") return(peak)
    peak1 <- as.data.frame(peak)
    if(!any(peak1$strand %in% c("+", "_")) && "feature_strand" %in% names(peak1)) peak1$strand <- peak1$feature_strand
    if(!any(peak1$strand %in% c("+", "_"))) {
        peak2 <- peak1
        peak1$strand <- "+"
        peak2$strand <- "-"
        peak1 <- rbind(peak1, peak2)
    }
    write.table(peak1[,c("name", "seqnames", "start", "end", "strand"), drop=FALSE], file=file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    return(invisible(peak))
} #}}}
my.homer <- function(bed.file, out.dir, run=FALSE, genome=NULL){ #{{{
    if(exists("GENOME") && is.null(genome)) genome <- GENOME
    #system("module load homer")
    homer.dir <- file.path(REF_DIR, "chipseq_pipeline_genome_data", paste0(genome, ".homer"))
    o.dir <- setwd(dirname(bed.file)); on.exit(setwd(o.dir))
    cmd <- paste0("findMotifsGenome.pl ", basename(bed.file),  " ", genome, "  ", out.dir,  " -size 200  -preparsedDir ", homer.dir)
    print(cmd)
    if(run) system(cmd)
} #}}}
my.homer.with.selfdefined.motifs <- function(bed.file, out.dir, run=FALSE, genome=NULL, motif.file="/scratch/gtac/analysis/chip_seq/bednarski_chipseq/homer/self_defined_ETS.motifs"){ #{{{
    if(exists("GENOME") && is.null(genome)) genome <- GENOME
    #system("module load homer")
    homer.dir <- file.path(REF_DIR, "chipseq_pipeline_genome_data", paste0(genome, ".homer"))
    o.dir <- setwd(dirname(bed.file)); on.exit(setwd(o.dir))
    cmd <- paste0("findMotifsGenome.pl ", basename(bed.file),  " ", genome, "  ", out.dir,  " -size 200  -preparsedDir ", homer.dir, " -mknown ",
                  motif.file, " -nomotif ") # search motifs from motif.file, and disable search for de novo motifs
    print(cmd)
    if(run) system(cmd)
} #}}}
find.peak.for.motif <- function(peak.txt, homer.out.dir, run=FALSE, genome=NULL){ #{{{
    if(exists("GENOME") && is.null(genome)) genome <- GENOME
    peak.txt0 <- normalizePath(peak.txt)
    peak.txt <- basename(peak.txt)
    o.dir <- setwd(homer.out.dir)
    on.exit(setwd(o.dir))
    bed <- gsub("txt$", "bed", peak.txt)
    fasta <- gsub("txt$", "fasta", peak.txt)
    sys.cmd("cut -f2-  ", peak.txt0, " > ", bed)
    fa <- get.ref.fa(genome)


    sys.cmd("bedtools getfasta -fi ", fa,
            " -bed ", bed, " -fo ", fasta)
    sys.cmd("cat ",
            paste(paste0("homerResults/motif", 1:10, ".motif"), collapse=' '),
            " > ", "top10.motif")
    cmd <- paste0("homer2 find -m top10.motif -i ", fasta, " > top10.motifs.and.peaks ")
    cat("cd", getwd(), "\n")
    cat(cmd, "\n")
    if(run) system(cmd)

    sys.cmd("cat ",
            paste(paste0("knownResults/known", 1:10, ".motif"), collapse=' '),
            " > ", "top10.known.motif")
    cmd <- paste0("homer2 find -m top10.known.motif -i ", fasta, " > top10.known.motifs.and.peaks ")
    cat(cmd, "\n")
    if(run) system(cmd)
    #cmd <- paste0("homer2 find -m /scratch/gtac/analysis/chip_seq/bednarski_chipseq/homer/self_defined_ETS.motfif -i ", fasta, " > self_defined_ETS.motifs.and.peaks ")
    #cat(cmd, "\n")
} #}}}
add.motif <- function(peak.name, motif.dir){ #{{{
    if("GRanges" %in% class(peak.name))
        peak.name <- paste0(seqnames(peak.name), ":", start(peak.name), "-", end(peak.name))
    motif <- rbind(read.delim(file.path("homer", motif.dir, "top10.known.motifs.and.peaks"), head=FALSE, as.is=TRUE),
                   read.delim(file.path("homer", motif.dir, "self_defined_ETS.motifs.and.peaks"), head=FALSE, as.is=TRUE))
    ret <- sapply(peak.name, function(x) paste(motif$V4[motif$V1 %in% x], collapse="; "))
    return(ret)
} #}}}
my.homer1liner <- function(peak, file, out.dir=file, rm.dup=TRUE, genome=NULL){ #{{{
    bed.file <- paste0(file, ".peak4homer.txt")
    write.peak4homer(peak, file=bed.file, rm.dup=rm.dup)
    if(!file.exists(out.dir)) dir.create(out.dir)
    out.dir <- normalizePath(out.dir)
    my.homer(bed.file, out.dir, run=TRUE, genome=genome)
    find.peak.for.motif(bed.file, out.dir, run=TRUE, genome=genome)
} #}}}
my.homer.example <- function(){ #{{{
    at('### example code for running homer
    write.peak4homer(shared.enhancer, file="slide1.c.shared.enhancer.peak4homer.txt", rm.dup=TRUE)
    my.homer("slide1.c.shared.enhancer.peak4homer.txt", "slide1.c.homer", run=TRUE, genome="hg19")
    find.peak.for.motif("slide1.c.shared.enhancer.peak4homer.txt", "slide1.c.homer", run=TRUE, genome="hg19")\n')
} #}}}
