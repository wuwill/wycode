find.file <- function(path, pattern, pattern2=NULL, vpattern=NULL){ #{{{
    cmd <- paste0("find ", path, " -name '", pattern, "'")
    if(!is.null(pattern2)) cmd <- paste0(cmd, " | grep '", pattern2, "'")
    if(!is.null(vpattern)) cmd <- paste0(cmd, " | grep -v '", vpattern, "'")
    return(system(cmd, TRUE))
} #}}}
get.new.chip <- function(path, type, pooled=TRUE, ...){ #{{{
    # type in bam, bw, peak
    ret <-
        if(pooled){ #{{{
            if(type %in% "peak") find.file(path, "optimal_peak.narrowPeak.gz") else
                if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="pooled") else 
                    find.file(path, paste0("*", type), ...)
        } else { # for replicates
            if(type %in% "peak") find.file(path, "*bfilt.narrowPeak.gz", pattern2="call-macs2/") else
                if(type %in% "bam") find.file(path, "*nodup.bam", pattern2="call-filter/", vpattern="inputs") else 
                    if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="call-macs2/", vpattern="inputs") else
                        find.file(path, paste0("*", type), ...)
        } #}}}
    return(ret)
} #}}}
get.new.atac <- function(path, type, pooled=TRUE, ...){ #{{{
    # type in bam, bw, peak
    ret <-
        if(pooled){ #{{{
            if(type %in% c("peak", "idr", "overlap")) grep(type, find.file(path, "optimal_peak.narrowPeak.gz", "execution/"), value=TRUE) else
                if(type %in% "bw") find.file(path, "*pooled.fc.signal.bigwig") else 
                    find.file(path, paste0("*", type), "execution/", ...)
        } else { # for replicates
            if(type %in% "peak") find.file(path, "*bfilt.narrowPeak.gz", pattern2="call-macs2/") else
                if(type %in% "bam") find.file(path, "*nodup.bam", pattern2="call-filter/", vpattern="inputs") else 
                    if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="call-macs2.*/shard", vpattern="inputs") else
                        find.file(path, paste0("*", type), ...)
        } #}}}
    return(ret)
} #}}}
get.prev.chip <- function(path, type, pooled=TRUE, ...){ #{{{
    # type in bam, bw, peak
    ret <-
        if(pooled){ #{{{
            if(type %in% "peak") find.file(path, "*filt.narrowPeak.gz", pattern2="optimal_set") else
                if(type %in% "bw") {
                    f <- find.file(path, "*fc.signal.bw")
                    if(length(f)>1) f <- grep("pooled_rep", f, value=TRUE)
                    f
                } else 
                    find.file(path, paste0("*", type), ...)
        } else { # for replicates
            if(type %in% "peak") find.file(path, "*filt.narrowPeak.gz", pattern2="macs2/rep") else
                if(type %in% "bam") find.file(path, "*.markdup.bam") else  #@ this does not work for now
                    if(type %in% "bw") find.file(path, "*fc.signal.bw", pattern2="macs2/rep", vpattern="inputs") else
                        find.file(path, paste0("*", type), ...)
        } #}}}
    return(ret)
} #}}}
get.prev.atac <- function(path, type, pooled=TRUE, ...){ #{{{
    # type in bam, bw, peak
    ret <-
        if(pooled){ #{{{
            if(type %in% "peak") find.file(path, "*narrowPeak.gz", pattern2="optimal") else
                if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="pooled_") else 
                    find.file(path, paste0("*", type), ...)
        } else { # for replicates
            if(type %in% "peak") find.file(path, "*filt.narrowPeak.gz", pattern2="macs2/rep") else
                if(type %in% "bam") find.file(path, "*nodup.bam", pattern2="align/rep") else 
                    if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="macs2/rep") else
                        find.file(path, paste0("*", type), ...)
        } #}}}
    return(ret)
} #}}}
get.chip.file <- function(group, type, ..., pattern=NULL){ #{{{
    i <- match(group, chip.dir$group)
    path <- file.path(chip.dir$top.path[i], chip.dir$dir[i])
    fun <- if(chip.dir$top.path[i] %in% new.pipe.list) get.new.chip else get.prev.chip
    ret <- fun(path, type, ...)
    if(!is.null(pattern)) grep(pattern, ret, value=TRUE) else return(ret)
} #}}}
get.atac.file <- function(group, type, ..., pattern=NULL){ #{{{
    i <- match(group, atac.dir$group)
    path <- file.path(atac.dir$top.path[i], atac.dir$dir[i])
    fun <- if(atac.dir$top.path[i] %in% new.pipe.list) get.new.atac else get.prev.atac
    ret <- fun(path, type, ...)
    if(!is.null(pattern)) grep(pattern, ret, value=TRUE) else return(ret)
} #}}}
get.file <- function(..., atac=TRUE){ #{{{
    if(atac) get.atac.file(...) else get.chip.file(...)
} #}}}
