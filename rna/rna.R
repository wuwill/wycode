library(Hmisc)
library(readxl)
library(heatmap3)
options(bitmapType="cairo")
get_colors <- function(groups, legend = FALSE) {#  {{{

#Custom Color Palettes#
	color.pals <- list(
		g1qualitative=c("#4477AA"),
		g2qualitative=c("#4477AA", "#CC6677"),
		g3qualitative=c("#4477AA", "#DDCC77", "#CC6677"),
		g4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
		g5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
		g6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
		g7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
		g8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499"),
		g9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
		g10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
		g11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
		g12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"),
		g13qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499", "#771122"),
		g14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
		g15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
		#g18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
		#g21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
	 )
  	 groups <- as.factor(groups)
  	 ngrps <- length(levels(groups))
	 group.col <- color.pals[[ngrps]]
    	 group.col <- rep(group.col, ngrps)
  	 color <- group.col[as.numeric(groups)]
  	 names(color) <- as.vector(groups)
	
	 if (legend == "TRUE") {
	   return(group.col)
	 } else {
  	   return(color)
	 }
} #}}}
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
get.diff <- function(group, project.dir=".", suffix=c('tsv', "xls")){ #{{{
    if(project.dir!=".") {
        o.dir <- setwd(project.dir); on.exit(setwd(o.dir))
    }
    library(readr)
    ret.file <- Sys.glob(paste0(group, "/*_Differential_Expression.", suffix))
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
plot_Volcano <- function(x, results_name, highlight_genes=NULL, alpha=0.6, i_significant=NULL, col_sig = "red", text_col = "black", text_size = 2.5, vline = NULL, hline = NULL, height = 1200, width = 1500, res = 300, ylim = NULL) {#  {{{
	require(ggplot2)
    library(ggrepel)
    x$log10p <- -log10(x$P.Value)
    if(is.null(i_significant)) {
        i_significant <- x$adj.P.Val <= 0.05
    }
    x$Sig <- "NS"
    if(is.list(i_significant)) {#{{{
        for(i in seq_along(i_significant)) {
            name_ <- names(i_significant)[i]
            x$Sig[i_significant[[i]]] <- name_
        }
        names(col_sig) <- names(i_significant)
        col_sig_ <- if("NS" %in% x$Sig) c(NS = "gray", col_sig) else col_sig
    } else if(is.numeric(i_significant) || is.logic(i_significant)) {
        x$Sig[i_significant] <- "Significant"
        col_sig_ <- c("NS" = "gray", Significant = col_sig)
    } else if(is.character(i_significant) || is.factor(i_significant)) {
        x$Sig <- i_significant
        col_sig_ <- col_sig
    }#}}}
    x$label <- ""
    if(!is.null(highlight_genes)) {
        i <- match(highlight_genes, x$external_gene_name)
        x$label[i] <- x$external_gene_name[i]
    }

	xlimit <- max(c(max(x$logFC), abs(min(x$logFC))))
    ylimit <- if(is.null(ylim)) max(-log10(x$P.Value)) else ylim
    png0 <- paste(results_name, "Volcano_plot.png", sep="_")
    png1 <- paste(results_name, "Volcano_plot.with_gene_labels.png", sep="_")
    # png3 <- paste(results_name, "Volcano_plot.with_gene_labels2.png", sep="_")
	p <- ggplot(x, aes(x = logFC, y = -log10(P.Value), colour = Sig, label = label)) + geom_point(size = 1.5, alpha = alpha) + xlim(-xlimit, xlimit) + scale_colour_manual(values = col_sig_) + xlab("Log 2 Fold-Change") + ylim(0, ylimit)
    if(length(col_sig_)>0) {
        for(i in 2:length(col_sig_))
            p <- p + geom_point(data = subset(x, Sig == names(col_sig_)[i]),
                                aes(x = logFC, y = -log10(P.Value), colour = Sig),
                                           size = 1.5, alpha = alpha)
    }
    if(!is.null(vline)) p <- p + geom_vline(xintercept = vline, colour = "gray", size = 1, alpha = 0.3)
    if(!is.null(hline)) p <- p + geom_hline(yintercept = vline, colour = "gray", size = 1, alpha = 0.3)

    png(png0, height=height, width=width, res = 300); print(p + theme_classic()); dev.off()

    if(!is.null(highlight_genes)){ #{{{
        p1 <-  p + geom_text_repel(box.padding = 0.5, max.overlaps = Inf,
                                   size = text_size,
                                   min.segment.length = 0,
                                   colour = text_col)


        # p1 <- geom_text_repel(data=x2, aes(x = logFC, y = log10p, label=external_gene_name), colour=text_col,
                              # size=2.5, min.segment.length = 0, segment.size = 0.8,
                              # box.padding = 0.5, max.overlaps = Inf)
        png(png1, height=height, width=width, res = 300); print(p1+theme_classic()); dev.off()
    } #}}}
} #}}}

plot_heatmap <- function(dat, sample.info, out, ann_gene = NULL, group_order = NULL, width = 4, height = 6, ...) {#{{{
    library(ComplexHeatmap)
    dat <- sweep(dat, 1, apply(dat, 1, mean, na.rm = TRUE))
    dat <- sweep(dat, 1, apply(dat, 1, sd, na.rm = TRUE), "/")
    ngenes <- nrow(dat)
    # genes <- ann$external_gene_name[match(rownames(dat), ann$Feature_ID)]
    genes <- rownames(dat)
    if(is.null(group_order)) group_order <- sort(unique(sample.info$group))

    # if(ngenes > 300) return()
    fs <- if(ngenes < 10) 12 else 
            round(10 / (sqrt(ngenes) - 3)) + 2

    pdf(paste0(out, ".heatmap.pdf"), width = width, height=height); on.exit(dev.off())

    if(!is.null(ann_gene) && any(ann_gene %in% genes)) {
        ann_gene <- ann_gene[ann_gene %in% genes]
        at <- match(ann_gene, genes)
        ha <- rowAnnotation(foo = anno_mark(at = at, labels = ann_gene))
    } else ha <- NULL
    show_row_names <- is.null(ann_gene) & !is.null(genes)
    # ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = month.name[1:10]))
    # Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha)
    print(Heatmap(as.matrix(dat), name = " ",
            column_split = factor(as.character(sample.info$group), levels = group_order),
            show_row_names = show_row_names,
            right_annotation = ha,
            row_names_gp = gpar(fontsize = fs),
            cluster_row_slices = FALSE, 
            cluster_column_slices = FALSE,
            ...))
}#}}}

get.gene_count <- function(key.file, ann.file="all.gene_counts.tsv", top.dir=".", count.type="RPKM"){#  {{{
    #key.file <- "" # go to each key and use "gene_counts.txt
    #ann.file <- "" # get annotation from this 
    key <- read.delim(key.file, head=FALSE, as.is=TRUE, sep="\t")
    names(key) <- c("Index", "index", "name")
    key$name <- paste0("sample.", key$name)

    anno <- read.delim(ann.file, head=TRUE, as.is=TRUE, se="\t")[1:7]

    ret <- sapply(1:nrow(key), function(i){
                  input <- read.delim(file.path(top.dir, key$index[i], "gene_counts.txt"), head=TRUE, as.is=TRUE)
                  input[match(anno[[1]], input[[1]]), count.type]
                    })
    colnames(ret) <- key$name
    return(cbind(anno, ret))
}#}}}

# pathway
get.fisher.p <- function(i, n1, n2, n){ #{{{
    a <- i
    b <- n2 - i
    c <- n1 - i
    d <- (n - n1) - (n2 -i)
    my.fisher <- function(a, b, c, d){
        dat <- matrix(c(a, b, c, d), nrow=2)
        fisher.test(dat, alternative="greater")$p.value
    }
    return(my.fisher(a, b, c, d))
} #}}}
get_msigdbr_dt <- function(species, pathway_types = c("GO_BP", "GO_MF", "KEGG")) {#{{{
    library(msigdbr)
    library(clusterProfiler)
    species_map <- c("mouse" = "Mus musculus",
                     "human" = "Homo sapiens",
                     "fly" = "Drosophila melanogaster",
                     "rat" =  "Rattus norvegicus")
    species <- plyr::mapvalues(species, names(species_map), species_map, warn_missing = FALSE)
    # msigdbr_show_species() # c("Bos taurus", "Caenorhabditis elegans", "Canis lupus familiaris", "Danio rerio", "Drosophila melanogaster", "Gallus gallus", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Saccharomyces cerevisiae", "Sus scrofa")

    find_category <- function(pathway_type) {#{{{
        if(pathway_type == "KEGG") return(c("C2", "CP:KEGG"))
        res <- stringr::str_split(pathway_type, "_")[[1]]
        if(res[1] == "GO") res[1] <- "C5"
        return(res)
    }#}}}

    res <- list()
    for(pathway_type in pathway_types) {#{{{
        category <- find_category(pathway_type)
        if(length(category) == 1) {#{{{
            category <- category[1]
            subcategory <- NULL
        } else {
            subcategory <- category[2]
            category <- category[1]
        }#}}}
        subcategory <- if(length(category))
        res[[pathway_type]] <- msigdbr(species = species,
                                       category = category,
                                       subcategory = subcategory)
    }#}}}
    
    if(length(pathway_types) == 1) return(res[[1]]) else return(res)
}#}}}
subset_enrichResult = function(x, i) { #  {{{
    rownames(x@result) = x@result$ID

    x@result = x@result[i, , drop = FALSE]

    x@geneSets = x@geneSets[x@result$ID]
    return(x)
} #}}}
clusterProfiler_enricher <- function(input, gene_col = "", gene_type = "ENTREZID", species = "mouse", file_prefix = "enricher") {# over representation analysis {{{
    # input: a data frame of DGE analysis results with these columns
    #       - gene
    library(ggplot2)
    library(dplyr)
    library(clusterProfiler)

    pathway_dt <- get_msigdbr_dt(species)
    mk_fn <- function(pathway_name="", ...) {#{{{
        res <- paste0(file_prefix, ".", pathway_name, "/", basename(file_prefix), ".", pathway_name, ...)
        if(!file.exists(dirname(res))) dir.create(dirname(res), recursive = TRUE)
        return(res)
    }#}}}
    # load orgdb {{{
    if(tolower(species) %in% c("mouse", "mus musculus")) {
        library(org.Mm.eg.db)
        orgdb <- org.Mm.eg.db
    }
    if(tolower(species) %in% c("human", "homo sapiens")) {
        library(org.Hs.eg.db)
        orgdb <- org.Hs.eg.db
    }
    if(tolower(species) %in% c("rat", "rattus norvegicus")) {
        library(org.Rn.eg.db)
        orgdb <- org.Rn.eg.db
    }
    if(tolower(species) %in% c("fly", "drosophila melanogaster")) {
        library(org.Dm.eg.db)
        orgdb <- org.Dm.eg.db
    }
    #}}}
    # convert input to ENTREZID vector {{{
    if(!is.null(gene_col) || ! gene_col %in% c(NA, ""))
        input <- input[[gene_col]]
    if(!tolower(gene_type) %in% c("", "entrezid", "entrez", "entrezgene"))
        input <- bitr(input, fromType = gene_type, toType = "ENTREZID", OrgDb = orgdb)
    #}}}
    test_ <- function(input, db, pathway_type) {#{{{
        db <- db[, c('gs_name', "entrez_gene")]
        res <- enricher(input, TERM2GENE = db, pvalueCutoff = 1, qvalueCutoff = 1)
        if(is.null(res)) return(NULL)
        res <- setReadable(res, orgdb, keyType = "ENTREZID")
        write_tsv(as.data.frame(res), mk_fn(pathway_type, "over_representation_test.tsv"))
        try_plot_top_fdr <- function(res, fdr=0.05) {#{{{
            if(is.null(fdr) || is.na(fdr)) {#{{{
                res <- if(nrow(res) > 20) subset_enrichResult(res, 1:20) else res
                name1 <- "top20_pathways"
                name2 <- "Top 20 pathways"
            } else {
                res <- subset_enrichResult(res, res@result$qvalue <= fdr)
                name1 <- paste0("pathway_significant_at_FDR", fdr)
                name2 <- paste0("FDR <= ", fdr)
            }#}}}
            n <- nrow(res)
            name <-
            if(n>1 && n <= 40) {
                p <- dotplot(res, showCategory = n, font.size = 9) + ggtitle(name2)
                ggsave(mk_fn(pathway_type, name1, ".dotplot.pdf"),
                       p,
                       width = 8,
                       height = max(8, round(0.8 * n) / 2))
            }
        }#}}}
        for(fdr in c(NA, 0.01, 0.05, 0.1, 0.2)) try_plot_top_fdr(res, fdr)
        return(res)
    }#}}}

    if(is.data.frame(pathway_dt[[1]])) {
        res <- list()
        for(pathway_type in names(pathway_dt)) 
            res[[pathway_type]] <- test_(input, pathway_dt[[pathway_type]], pathway_type)
    } else {
        res <- test_(input, pathway_dt, pathway_type)
    }

    return(res)
}#}}}
