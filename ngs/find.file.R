find.file <- function(path, pattern, pattern2=NULL, vpattern=NULL){ #{{{
    cmd <- paste0("find ", path, " -name '", pattern, "'")
    if(!is.null(pattern2)) cmd <- paste0(cmd, " | grep '", pattern2, "'")
    if(!is.null(vpattern)) cmd <- paste0(cmd, " | grep -v '", vpattern, "'")
    return(system(cmd, TRUE))
} #}}}
get.chip.file <- function(path, type, pooled=TRUE){ #{{{
    # type in bam, bw, peak
    ret <-
        if(pooled){ #{{{
            if(type %in% "peak") find.file(path, "optimal_peak.narrowPeak.gz") else
                if(type %in% "bw") find.file(path, "*pooled.fc.signal.bigwig") else
                    ""
        } else { # for replicates
            if(type %in% "peak") find.file(path, "*bfilt.narrowPeak.gz", pattern2="call-macs2/") else
                if(type %in% "bam") find.file(path, "*nodup.bam", pattern2="call-filter/", vpattern="inputs") else
                    if(type %in% "bw") find.file(path, "*fc.signal.bigwig", pattern2="call-macs2/", vpattern="inputs") else
                        ""
        } #}}}
    return(ret)
} #}}}
