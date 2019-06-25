library(Hmisc)
library(readxl)
library(heatmap3)
options(bitmapType="cairo")
venn4de <- function(groups, file, project.dir="."){ # read signficant DE resutls and make Venn Diagram {{{
    o.dir <- setwd(project.dir); on.exit(setwd(o.dir))
    get.sig <- function(group){ #{{{
        library(readr)
        ret.file <- Sys.glob(file.path(group, "*_low_FDR_Significant_Differential_Expression_Subset.xls"))
        ret <- read_tsv(ret.file)
        return(ret[[1]])
    } #}}}
    ret <- lapply(groups, get.sig)
    names(ret) <- groups
    library(VennDiagram)
    venn.diagram(ret, euler.d = TRUE, scaled=TRUE, filename=file,
                 alpha = 0.50,
                 col = "transparent",
                 fill = c("cornflowerblue", "green", "yellow")#, "darkorchid1")
                 )
    venn.diagram(ret, euler.d = TRUE, scaled=TRUE,filename=gsub(".tiff$", ".no_text.tiff", file), 
                 alpha = 0.50,
                 col = "transparent",
                 fill = c("cornflowerblue", "green", "yellow")#, "darkorchid1")
                 )
    venn.diagram(ret, euler.d = FALSE, scaled=FALSE, filename=gsub(".tiff$", ".fixed_size.tiff", file),
                 alpha = 0.50,
                 col = "transparent",
                 fill = c("cornflowerblue", "green", "yellow")#, "darkorchid1")
                 )
    venn.diagram(ret, euler.d = FALSE, scaled=FALSE, filename=gsub(".tiff$", ".fixed_size.no_text.tiff", file), 
                 alpha = 0.50,
                 col = "transparent",
                 category.names=rep("", 3),
                 fill = c("cornflowerblue", "green", "yellow")#, "darkorchid1")
                 )
    return(invisible(ret))
} #}}}
read.rna.result <- function(dex.file=NULL){ # get dat / sample.info / ann{{{
    # expression data
    if(!is.null(dex.file)) load(dex.file)
    dat <- d$E
    # sample.info
    sample.info <- d$targets
    names(sample.info)[3] <- "Name"
    sample.info$Name <- gsub("^sample.", "", sample.info$Name)
    genes <- rownames(d)

    i.ann <- match(genes, annotations$Feature_ID)
    ann <- annotations[i.ann,]
    ann$Feature_ID <- genes

    colnames(dat) <- sample.info$Name
    dat <<- dat
    sample.info <<- sample.info
    ann <<- ann
    fit <<- fit

    return(invisible())
} #}}}
plotHeatmap <- function(y, genes, heatmap_name, groups, main=NULL, ...) {#  {{{
	    y <- y[!duplicated(genes),,drop=FALSE]
	    rownames(y) <- genes
	    plot_height <- min(c(max(c(1500, (nrow(y)*10))), 20000))
	    plot_width <- min(c(max(c(1500, (nrow(y)*15))), 20000))
	    
	    if (nrow(y) <= 30) {
	       textsize <- 1
	    }
	    else if ((nrow(y) > 30) & (nrow(y) <= 1500)) {
	    	textsize <- 0.8
	    } else {
	      	textsize <- 0.1
	    }
	    
        if(is.null(main)) main <- basename(heatmap_name)
        heatmap_name <- file.path(dirname(heatmap_name), gsub("/", ":", basename(heatmap_name)))
	    png(paste0(heatmap_name, ".png"), height = plot_height, width = plot_width)#, type = "cairo")
        heatmap3(y, balanceColor = T, ColSideAnn=as.data.frame(groups), ColSideFun=function(x) showAnn(x), ColSideWidth = 0.8, main = main, xlab = "Samples", ylab = "Genes", margins = c(15,10), cexRow = textsize, ...)
		#heatmap3(y, balanceColor = T, main = heatmap_name, xlab = "Samples", ylab = "Genes", margins = c(15,10), cexRow = textsize, ...)
	    dev.off()
}#}}}
my.plot.heatmap.pair <- function(group.pair, project.dir=".", genes=NULL, suffix=NULL, gene.id="external_gene_name", Colv=NULL, ...){#  {{{
    o.dir <- setwd(project.dir); on.exit(setwd(o.dir))
    if(length(group.pair)==1) group.pair <- strsplit(group.pair, "__vs__")[[1]]
    group <- paste(group.pair, collapse="__vs__")
    get.sig <- function(group){ #{{{
        library(readr)
        ret.file <- Sys.glob(file.path(group, "*_low_FDR_Significant_Differential_Expression_Subset.xls"))
        ret <- read_tsv(ret.file)
        return(setdiff(ret[["external_gene_name"]], c(NA, "")))
    } #}}}
    if(is.null(genes)) genes <- get.sig(group) else { #{{{
        genes <- setdiff(ann$external_gene_name[match(genes, ann[[gene.id]])], NA)
    } #}}}
    if(length(genes)<2) {
        cat(suffix, "- Skipped because number of genes with external_gene_name < 2\n")
        return()
    }
    d <- dat[match(genes, ann$external_gene_name), sample.info$group %in% group.pair, drop=FALSE]
    samples <- colnames(d)
    groups <- sample.info$group[match(samples, sample.info$Name)]
    groups <- as.character(groups)
    i <- order(groups)
    file.name <- if(is.null(suffix)) group else file.path(dirname(suffix), paste0(group, ".", basename(suffix)))
    #plotHeatmap(d, genes, paste0(file.name, ".heatmap"), groups)
    plotHeatmap(d[,i], genes, paste0(file.name, ".byGroup.heatmap"), groups[i], Colv=Colv)
}#}}}
get.diff <- function(group, project.dir="."){ #{{{
    if(project.dir!=".") {
        o.dir <- setwd(project.dir); on.exit(o.dir)
    }
    library(readr)
    ret.file <- Sys.glob(file.path(group, "*_Differential_Expression.xls"))
    ret <- read_tsv(ret.file)
    return(as.data.frame(ret))
} #}}}
plot_MA <- function(DIFF, results_name, highlight_genes=NULL, alpha=0.6) { #{{{
    png0 <- paste(results_name, "MA_plot.png", sep="_")
    png1 <- paste(results_name, "MA_plot.no_labels.png", sep="_")
    png2 <- paste(results_name, "MA_plot.with_gene_labels1.png", sep="_")
    png3 <- paste(results_name, "MA_plot.with_gene_labels2.png", sep="_")
    p <- ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value, label=external_gene_name)) + geom_point(size = 4, alpha = 0.5) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot") #+ 
    #geom_text(aes(label=ifelse(external_gene_name %in% highlight_genes,as.character(external_gene_name),'')),hjust=0.5,vjust=0)
    if(!is.null(highlight_genes)){ #{{{
        i <- match(highlight_genes, DIFF$external_gene_name)
        diff2 <- DIFF[i,,drop=FALSE]
        diff1 <- DIFF[-i,,drop=FALSE]
        p <- ggplot(diff1, aes(x = AveExpr, y = logFC, colour = P.Value, label=external_gene_name)) + geom_point(size = 4, alpha = 0.5) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot") #+ 
        p1 <- geom_point(data=diff2, aes(x = AveExpr, y = logFC),colour="black", size=4)
        p1.5 <- geom_point(data=diff2, aes(x = AveExpr, y = logFC),colour="black", size=4, alpha=alpha)
        p2 <- geom_text(data=diff2, aes(x = AveExpr, y = logFC, label=external_gene_name), colour="black", hjust=-0.33, vjust=0.5, size=4)
        p3 <- geom_text(data=diff2, aes(x = AveExpr, y = logFC, label=external_gene_name), colour="black", hjust=0.5, vjust=-1, size=4)
        png(png1, height=600, width=600); print(p+p1); dev.off()
        png(png2, height=600, width=600); print(p+p1.5+p2); dev.off()
        png(png3, height=600, width=600); print(p+p1.5+p3); dev.off()
    } else {
        png(png0, height=600, width=600); print(p); dev.off()
    } #}}}
} #}}}
plot_Volcano <- function(x, results_name, highlight_genes=NULL, alpha=0.6) {#  {{{
	require(ggplot2)
    x$log10p <- -log10(x$P.Value)
    DIFF <- x
	xlimit <- max(c(max(x$logFC), abs(min(x$logFC))))
    png0 <- paste(results_name, "Volcano_plot.png", sep="_")
    png1 <- paste(results_name, "Volcano_plot.no_label.png", sep="_")
    png2 <- paste(results_name, "Volcano_plot.with_gene_labels1.png", sep="_")
    png3 <- paste(results_name, "Volcano_plot.with_gene_labels2.png", sep="_")
	p <- ggplot(x, aes(x = logFC, y = -log10(P.Value), colour = adj.P.Val <= 0.05)) + geom_point(size = 1.5, alpha = 0.7) + xlim(-xlimit, xlimit) + geom_vline(xintercept = 2, colour = "blue", size = 1, alpha = 0.3) + geom_vline(xintercept = -2, colour = "blue", size = 1, alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) + xlab("Log 2 Fold-Change")
    if(!is.null(highlight_genes)){ #{{{
        i <- match(highlight_genes, x$external_gene_name)
        diff2 <- DIFF[i,,drop=FALSE]
        diff1 <- DIFF[-i,,drop=FALSE]
        p <- ggplot(diff1, aes(x = logFC, y = log10p, colour = adj.P.Val <= 0.05)) + geom_point(size = 1.5, alpha = 0.7) + xlim(-xlimit, xlimit) + geom_vline(xintercept = 2, colour = "blue", size = 1, alpha = 0.3) + geom_vline(xintercept = -2, colour = "blue", size = 1, alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) + xlab("Log 2 Fold-Change")
        p1 <- geom_point(data=diff2, aes(x = logFC, y = log10p),colour="darkblue", size=1.5)
        p1.5 <- geom_point(data=diff2, aes(x = logFC, y = log10p),colour="darkblue", size=1.5, alpha=alpha)
        p2 <- geom_text(data=diff2, aes(x = logFC, y = log10p, label=external_gene_name), colour="darkblue", hjust=-0.33, vjust=0.5, size=4)
        p3 <- geom_text(data=diff2, aes(x = logFC, y = log10p, label=external_gene_name), colour="darkblue", hjust=0.5, vjust=-1, size=4)
        png(png1, height=600, width=600); print(p+p1); dev.off()
        png(png2, height=600, width=600); print(p+p1.5+p2); dev.off()
        png(png3, height=600, width=600); print(p+p1.5+p3); dev.off()
    } else {
        png(png0, height=600, width=600); print(p); dev.off()
    } #}}}
} #}}}
