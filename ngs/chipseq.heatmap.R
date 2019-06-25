library(EnrichedHeatmap)
read.bigwig2rle <- function(bigwig, peaks=NULL, bp=500, aggregate=TRUE, group=gsub(".bigwig", "", bigwig), fun=function(x) x, ...){ #{{{
    if(is.null(peaks)) signal <- import(bigwig, format="BigWig", as="Rle")
    if(!is.null(peaks)){
        if(!is.null(bp)) peaks <- resize(peaks, width=bp*2+1, fix="center")
        signal <- import(bigwig, selection=BigWigSelection(peaks), format="BigWig", as="Rle")
        if(!is.null(bp)){ #{{{
            pad <- Rle(0, max(end(peaks)))
            pad.left <- Rle(0, bp)
            for(chr in seqlevels(signal)){ #{{{
                signal[[chr]] <- c(pad.left, signal[[chr]], pad)
            } #}}}
            signal <- signal[GenomicRanges::shift(peaks, bp)]
        } else {
            signal <- signal[peaks]
        } #}}}
    }
    names(signal) <- group
    signal <- sapply(signal, function(x) fun(as.vector(x)))
    if(aggregate){ #{{{
        return(apply(signal, 1, mean))
    } #}}}
    return(signal)
} #}}}
my.generate.split.plots <- function(peaks, heat.mats, file="", signal.convert=function(x) x, split.by=peaks$group, peak.names, ...){ #{{{
    if(file != ""){ #{{{
        pdf(file)
        on.exit(dev.off())
    } #}}}
    heat.mats <- lapply(heat.mats, signal.convert)
    peak.names <- names(heat.mats)
    ugroups <- if(is.null(peak.names) || is.null(split.by)) NULL else c("both", peak.names)
    col_fun = colorRamp2(quantile(unlist(heat.mats), c(0, 0.95)), c("white", "red"))
    #get.col_fun <- function(x) if(grepl("k4$", x)) col_fun_k4 else if(grepl("k27$", x)) col_fun_k27 else col_fun
    get.col_fun <- function(...) col_fun
    #row_order <- do.call(c, by(-i3.peaks$Fold, split, order))
    #h1 <- apply(heat.mats[[1]], 1, max)
    #h2 <- apply(heat.mats[[2]], 1, max)
    #d <- h1 - h2
    #row_order <- do.call(c, lapply(ugroups,
                                   #function(xi){
                                       #i<- which(split %in% xi)
                                       #i[order(d[i])]
                                   #}))
    #row_order <- order(-i3.peaks$pos$Fold)
    if(!is.null(ugroups)){ #{{{
        col <- c("magenta", "red", "blue")
        names(col) <- ugroups
        split <- factor(peaks$group, levels=ugroups)
        hm_split <- Heatmap(split, show_row_names = FALSE, name = "group",
                            col = col,
                            show_column_names = FALSE, width = unit(2, "mm"))
    } #}}}
    title_cex <- gpar(cex=0.6)
    i.sel <- 1
    heatmap1 <- EnrichedHeatmap(heat.mats[[i.sel]],
                                top_annotation=NULL,
                                col=col_fun,
                                #row_order=row_order,
                                name="",
                                column_title=names(heat.mats)[i.sel],
                                column_title_gp=title_cex)
    heatmaps <- c(heatmap1, lapply(names(heat.mats)[-i.sel], function(x)
                                   EnrichedHeatmap(heat.mats[[x]],
                                                   top_annotation=NULL,
                                                   col=get.col_fun(x),
                                                   #row_order=row_order,
                                                   column_title=x,
                                                   show_heatmap_legend=FALSE,
                                                   column_title_gp=title_cex)))[order(c(i.sel, seq_along(heat.mats)[-i.sel]))]
    names(heatmaps) <- names(heat.mats)

    for(i in seq_along(heatmaps)){ #{{{
        p <- if(i==1) heatmaps[[1]] else p + heatmaps[[i]]
    } #}}}

    if(is.null(ugroups)){ #{{{
        ht_list = draw(p, cluster_rows = TRUE,
                       #row_order = row_order,
                       show_row_dend = FALSE,
                       heatmap_legend_side = "bottom", gap = unit(2, "mm"))

    } else{
        ht_list = draw(hm_split+p, cluster_rows = TRUE,
                       #row_order = row_order,
                       show_row_dend = FALSE,
                       split = split, heatmap_legend_side = "bottom", gap = unit(2, "mm"))

    }#}}}

} #}}}
mk.heatmap <- function(combined.peaks, bw.files, bam.files=NULL, out.suffix="", bp=1000, peak.names=c("peak1", "peak2"), ...){ #{{{
    heat.mat.overlap <- lapply(bw.files, get.heat.mat, peaks=combined.peaks, bp=bp/2)
    ov.bed <- paste0(out.suffix, ".peak.combined.bed")
    
    if(!is.null(bam.files)){ #{{{
        library(CoverageView)
        export(combined.peaks, ov.bed)
        coverages <- lapply(bam.files, CoverageBamFile)
        cov.matrix.ret <- lapply(coverages, function(x) cov.matrix(x, coordfile=ov.bed, no_windows=100, offset=0, normalization="rpm"))
        for(i in seq_along(heat.mat.overlap)){ #{{{
            heat.mat.overlap[[i]][] <- t(cov.matrix.ret[[i]])
        } #}}}
    } #}}}
    #names(heat.mat.overlap) <- peak.names
    names(heat.mat.overlap) <- gsub("_histone$", "", sapply(strsplit(bw.files, "/"), grep, pattern="histone$", value=TRUE))
    
    my.generate.split.plots(combined.peaks, heat.mat.overlap, file=gsub("bed$", "heatmap.pdf", ov.bed), signal.convert=function(x) x, peak.names=peak.names, ...)
    return(list(peak.overlap=combined.peaks, heat.mat=heat.mat.overlap))
} #}}}
my.mk.signal.histgram <- function(aggregated.signal, file="", i.sel=1:ncol(aggregated.signal), signal.convert=function(x) log2(x+1)){ #{{{
    if(file != ""){ #{{{
        pdf(file)
        on.exit(dev.off())
    } #}}}
    cls10 <- alpha(brewer.pal(10, 'Set1'), 0.6)
    aggregated.signal <- signal.convert(aggregated.signal[,i.sel,drop=FALSE])
    nbp <- (nrow(aggregated.signal) - 1)/2
    n <- ncol(aggregated.signal)
    if(n<1) return()
    col <- cls10[1:n]
    matplot(x=-nbp:nbp, aggregated.signal, xlab='', ylab='', type='l', col=col, lty=1, lwd=3, main="")
    legend.items <- colnames(aggregated.signal)
    if(is.null(legend.items)) legend.items <- 1:n
    legend('topright', legend=legend.items, col=col, lty=1, lwd=3, cex=0.8, border=0)
    matplot(x=-nbp:nbp, aggregated.signal, xlab='', ylab='', type='l', col=col, lty=1, lwd=3, main="")
} #}}}
