#####################################
##Load Options##
#####################################
suppressPackageStartupMessages(library(optparse))
options(bitmapType="cairo")

option_list <- list(#  {{{
make_option(c("-k", "--samplekey"), action="store", help = "Samplekey pair file with no header information.
For microarray data the first column is the filename, second column is the samplename, path of the file or - for current directory, and the remaining columns specific to the statistical model employed.
For RNAseq data the samplekey pair file is the same as the samplekey with an extra column for the name of the analytical directory or - for current analytical directory, and the remaining columns specific to the statistical model employed.
"),
make_option(c("-p", "--platform"), action = "store", help = "Microarray platform to be analyzed:
agilent - agilent data exported from Feature Extraction software.

illumina - illumina .idat files

affy - affymetrix .cel files

rnaseq - featureCount count files

NOTE:  The default platform is agilent."),
make_option(c("-c", "--chip_type"), action = "store", help = "Platform specific chip-type to be analyzed:
Human:
efg_agilent_sureprint_g3_ge_8x60k
efg_agilent_sureprint_g3_ge_8x60k_v2
efg_agilent_sureprint_g3_ge_8x60k_v3
efg_agilent_4x44k_v1
efg_agilent_4x44k_v2
affy_hc_g110
affy_hg_focus
affy_hg_u133_plus_2
affy_hta_2_0
affy_huex_1_0_st_v2
affy_hugene_1_0_st_v1
affy_hugene_2_0_st_v1
illumina_humanht_12_v4
illumina_humanwg_6_v3
NULL - Must provide custom gene annotations

Mouse:
efg_agilent_sureprint_g3_ge_8x60k
efg_agilent_wholegenome_4x44k_v1
efg_agilent_wholegenome_4x44k_v2
affy_moex_1_0_st_v1
affy_mogene_1_0_st_v1
affy_mouse430_2
affy_mouse430a_2
illumina_mousewg_6_v2
illumina_mouseref_8_v2
NULL -Must provide custom transcript annotations

RNA-seq:
voom
"),
make_option(c("-f", "--feature"), action="store", default = "gene", help = "Feature level to be analyzed and annotated:  gene or transcript
"),
make_option(c("-s", "--species"), action="store", help = "Species: mouse, human, rat, zebrafish, tomato, arabidopsis, malaria, cow, chicken, crypto, fly, worm, yeast, and NULL.  If NULL, there will be no KEGG or GO analysis and no annotations provided with the results.
"),
make_option(c("-e", "--estLatentFactors"), action="store", default = "FALSE", help = "Estimate the number of latent factors that contribute to uknown sources of variation and add them to the design model
"),
make_option(c("-b", "--paired"), action="store", default = "FALSE", help = "Samples consist of mostly paired subjects: TRUE or FALSE
For paired samples of unequivalent rank where not every sample has a corresponding matched sample for each condition, add an additional column to the pair file that contains the matched subject ID's and set this option to TRUE.  If the paired samples are of equivalent rank where every condition has a corresponding matched sample, then use the additive model and set this option to FALSE.  If estLatentFactors with SVA is set to TRUE, then paired should always be set to FALSE. The default is FALSE.
"),
make_option(c("-y", "--miRNA"), action="store", default = FALSE, help = "miRNA only analysis
"),
make_option(c("-q", "--transgenic"), action="store", default = FALSE, help = "Transgenic only analysis
"),
make_option(c("-x", "--custom_gene_annotations"), action="store", default = NULL, help = "Custom gene annotations
"),
make_option(c("-z", "--custom_transcript_annotations"), action="store", default = NULL, help = "Custom transcript annotations
"),
make_option(c("-m", "--model"), action = "store", default = "standard", help = "Statistical model to be used: standard, additive, time-series, interaction, or fit

standard - generalized-linear model to test for differences between two or more conditions with a single common control group.  Pair file should have columns for index_name, index, samplename, path, and group.  The control group should be labeled control.

additive - generalized-linear model to account for categorical factors or random-effects.  This is appropriate for controlling for batch effects, perfectly paired samples across all conditions, or accounting for gender, etc.  Pair file should have columns for file, samplename, path, group to be tested for, and an extra column for the factor to be accounted for.  Controls should be labeled control in the group column.  The controls in the factor column will be picked based on alphanumeric order.

interaction - generalized-linear model to account for factors that interact with treatment groups.  Pair file should have columns for index_name, index, samplename, path, treatment, and an extra column for the interaction factor.  Controls should be labeled control.

time-series - generalized-linear model to account for treatment vs placebo over time (a specialized interaction model).  Pair file should have columns for index_name, index, samplename, path, treatment, and an extra column for the time.  Controls should be labeled control.

fit - load data and fit model but contrasts are to be done manually with makeContrasts to create the contrast matrix for differential expression.

fit-additive - load data and fit model but contrasts are to be done manually with makeContrasts to create the contrast matrix for differential expression.

NOTE: To perform contrasts between coefficients, run the appropriate model to completion, open R, load the saved image, check the design model coefficients with as.data.frame(colnames(fit$coefficients)) and create a contrast.R Rscript with the limma, gage, pathview, ggplot2, and heatmap3 libraries, load the DEX_analysis.RData image, and create contrast matrices with the makeContrasts function.  Then execute the GLMmatrix function.  The GLMmatrix function will automatically generate Venn diagrams if the number of contrasts per contrast matrix is 5 or less.  The GLMcoefficients function can be used instead to test single coefficients of the design model or report changes across multiple coefficients, in which case the makeContrasts function is not required.
"))

opt <- parse_args(OptionParser(usage = "%prog [options]", option_list=option_list))

#####################################
##Load libraries and Option Parsing##
#####################################
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(gage))
suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(heatmap3))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(rgl))

r3dDefaults$windowRect = c(0,0,800,1000)

if (is.null(opt$samplekey)) {#  {{{
   cat("Missing samplekey pair file....try again\n\n")
   stop()
} else {
  samplekey <- opt$samplekey
}#}}}

if (opt$platform != "agilent" & opt$platform != "illumina" & opt$platform != "affy" & opt$platform != "rnaseq") {#  {{{
	cat("Platform must be agilent, illumina, affymetrix, or rnaseq......try again\n\n")
	stop()
} else {
       platform <- opt$platform
}#}}}

if (opt$chip_type != "efg_agilent_sureprint_g3_ge_8x60k"
& opt$chip_type != "efg_agilent_sureprint_g3_ge_8x60k_v2"
& opt$chip_type != "efg_agilent_sureprint_g3_ge_8x60k_v3"
& opt$chip_type != "efg_agilent_4x44k_v1"
& opt$chip_type != "efg_agilent_4x44k_v2"
& opt$chip_type != "affy_hc_g110"
& opt$chip_type != "affy_hg_focus"
& opt$chip_type != "affy_hg_u133_plus_2"
& opt$chip_type != "affy_hta_2_0"
& opt$chip_type != "affy_huex_1_0_st_v2"
& opt$chip_type != "affy_hugene_1_0_st_v1"
& opt$chip_type != "affy_hugene_2_0_st_v1"
& opt$chip_type != "illumina_humanht_12_v4"
& opt$chip_type != "illumina_humanwg_6_v3"
& opt$chip_type != "efg_agilent_wholegenome_4x44k_v1"
& opt$chip_type != "efg_agilent_wholegenome_4x44k_v2"
& opt$chip_type != "affy_moex_1_0_st_v1"
& opt$chip_type != "affy_mogene_1_0_st_v1"
& opt$chip_type != "affy_mouse430_2"
& opt$chip_type != "affy_mouse430a_2"
& opt$chip_type != "illumina_mousewg_6_v2"
& opt$chip_type != "illumina_mouseref_8_v2"
& opt$chip_type != "NULL"
& opt$chip_type != "voom") {
  cat("Chip type not supported..........try again\n\n")
  stop()
} else {
  chip_type <- opt$chip_type
}

if (opt$chip_type == "NULL" & opt$feature == "gene" & is.null(opt$custom_gene_annotations)) {
	cat ("Chip type is NULL, but custom gene annotations not provided......try again\n\n")
	stop()
}

if (opt$chip_type == "NULL" & opt$feature == "transcript" & is.null(opt$custom_transcript_annotations)) {
	cat ("Chip type is NULL, but custom transcript annotations not provided......try again\n\n")
	stop()
}

if (platform == "illumina") {
   if (chip_type == "illumina_humanht_12_v4") {
      bgx <- "/srv/seq/analysis1/microarray/illumina/HumanHT-12_V4_0_R2_15002873_B.bgx"
   } else if (chip_type == "illumina_mousewg_6_v2") {
      bgx <- "/srv/seq/analysis1/microarray/illumina/MouseWG-6_V2_0_R3_11278593_A.bgx"
   } else if (chip_type == "illumina_mouseref_8_v2") {
      bgx <- "/srv/seq/analysis1/microarray/illumina/MouseRef-8_V2_0_R3_11278551_A.bgx"
   } else if (chip_type == "illumina_humanwg_6_v3") {
      bgx <- "/srv/seq/analysis1/microarray/illumina/HumanWG-6_V3_0_R3_11282955_A.bgx"
   }
}

if (opt$species == "NULL") {
   go_species <- NULL
   kegg_species <- NULL
   biomart <- NULL
   reference_database <- NULL
   go_db <- NULL
   kegg_db <- NULL
}

cat ("Reading in options and databases...\n")

if (opt$species == "human") {
	go_species <- "human"
	kegg_species <- "human"
	host <- "www.ensembl.org"
	biomart <- "ENSEMBL_MART_ENSEMBL"
	reference_database <- "hsapiens_gene_ensembl"
	go_db <- go.gsets(species = go_species, id.type = "entrez")
	kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "mouse") {
	go_species <- "mouse"
	kegg_species <- "mouse"
	#host <- "www.ensembl.org"
    host <- "http://apr2018.archive.ensembl.org"
	biomart <- "ENSEMBL_MART_ENSEMBL"
	reference_database <- "mmusculus_gene_ensembl"
	go_db <- go.gsets(species = go_species, id.type = "entrez")
	kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "fly") {
	go_species <- "fly"
	kegg_species <- "dme"
	host <- "www.ensembl.org"
	biomart <- "ENSEMBL_MART_ENSEMBL"
	reference_database <- "dmelanogaster_gene_ensembl"
	go_db <- go.gsets(species = go_species, id.type = "entrez")
	kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "worm") {
   go_species <- "worm"
   kegg_species <- "cel"
   host <- "www.ensembl.org"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "celegans_gene_ensembl"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "zebrafish") {
   go_species <- "zebrafish"
   kegg_species <- "dre"
   host <- "www.ensembl.org"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "drerio_gene_ensembl"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "yeast") {
   go_species <- "yeast"
   kegg_species <= "sce"
   host <- "www.ensembl.org"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "scerevisiae_gene_ensembl"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "chicken") {
   go_species <- "chicken"
   kegg_species <- "gga"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "ggallus_gene_ensembl"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "rat") {
   go_species <- "rat"
   kegg_species <- "rno"
   host <- "www.ensembl.org"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "rnorvegicus_gene_ensembl"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "cow") {
   go_species <- "cow"
   kegg_species <- "bta"
   biomart <- "ENSEMBL_MART_ENSEMBL"
   reference_database <- "btaurus_gene_ensembl"
   go_db <- NULL  ##go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "arabidopsis") {
   go_species <- "arabidopsis"
   kegg_species <- "ath"
   host <- "plants.ensembl.org"
   biomart <- "plants_mart"
   reference_database <- "athaliana_eg_gene"
   go_db <- go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "tomato") {
   go_species <- NULL
   kegg_species <- "sly"
   host <- "plants.ensembl.org"
   biomart <- "ENSEMBL_MART_PLANT"
   reference_database <- "slycopersicum_eg_gene"
   go_db <- NULL  ###go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "crypto") {
   go_species <- NULL
   kegg_species <- "cne"
   host <- "plants.ensembl.org"
   biomart <- "fungi_mart_22"
   reference_database <- "cneoformans_eg_gene"
   go_db <- NULL  ####go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "entrez")
}

if (opt$species == "b_terrestris") {
   go_species <- NULL
   kegg_species <- "bter"
   host <- NULL
   biomart <- NULL
   reference_database <- NULL
   go_db <- NULL  ####go.gsets(species = go_species, id.type = "entrez")
   kegg_db <- kegg.gsets(species = kegg_species, id.type = "kegg")
}

if (opt$model != "additive" & opt$model != "time-series" & opt$model != "standard" & opt$model != "interaction" & opt$model != "fit" & opt$model != "fit-additive") {
	cat("Model must be standard, additive, interaction, time-series, or fit......try again\n\n")
	stop()
} else {
	model <- opt$model
}
#}}}

#####################################
##Functions##
#####################################
trim <- function (x) {gsub("^\\s+|\\s+$", "", x)}
fixnames <- function(x) { #  {{{
	 z <- gsub("\\s+", ".", perl = T, x)
	 gsub("-","_", perl = T, z)
} #}}}
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
plot_PCA <- function(x, name  = "PCA_plot", groups = NULL, all = FALSE, interactive = FALSE) { #  {{{
	 require(rgl)
	 require(RColorBrewer)
	 r3dDefaults$windowRect = c(0,0, 1200, 1600)
	 y <- x[,c(grep("sample", colnames(x)))]
	 colnames(y) <- gsub("sample.","", colnames(y))
	 title = "PCA Plot of all Features"
	 
	 if (all == FALSE) {
	    #Select for rank-variant genes with Std.Dev >= 1 to maximize clustering
	    row.sd <- apply(y,1,sd)
	    y <- y[row.sd >= 1,]
	    #row.var <- apply(y,1,var)
	    #y <- y[row.var >= 1,]
	    title = "PCA Plot of Variant Features with Std.Dev >= 1"
	 }
	 
	 p.comp <- prcomp(y)
	 p.comp.summary <- summary(p.comp)
	 p.comp.vectors <- as.data.frame(p.comp$rotation)
	 
	 plot3d(x = p.comp.vectors$PC1, xlab = paste("PC1:", p.comp.summary$importance[2,1]), y = p.comp.vectors$PC2, ylab = paste("PC2:", p.comp.summary$importance[2,2]), z = p.comp.vectors$PC3, zlab = paste("PC3:", p.comp.summary$importance[2,3]), type = "s", col =  get_colors(groups),  main = title)
	 legend3d("topright", levels(as.factor(groups)), col  = get_colors(groups, legend = TRUE), pch = 16, cex = 0.8, inset = c(0.02))
	 
	 if (interactive == TRUE) {
	    identify3d(x = p.comp.vectors$PC1,  y = p.comp.vectors$PC2, z = p.comp.vectors$PC3, labels = rownames(p.comp.vectors), buttons = c("right", "middle"))
	 } else {
	    text3d(x = p.comp.vectors$PC1,  y = p.comp.vectors$PC2, z = p.comp.vectors$PC3, texts = rownames(p.comp.vectors), adj = c(-0.5,1))
	    out.file <- writeWebGL(dir = "PCA", filename = paste(name,".html", sep = ""), width = 1200, height = 1600)
	    rgl.quit()
	}
} #}}}
plot_Volcano <- function(x, name = NULL) { #  {{{
	require(ggplot2)
	xlimit <- max(c(max(x$logFC), abs(min(x$logFC))))
	png(file = ifelse(is.null(name), "Volcano_plot.png", paste(name, "Volcano_plot.png", sep = "_")), height=600, width=600)
	print(ggplot(x, aes(x = logFC, y = -log10(P.Value), colour = adj.P.Val <= 0.05)) + geom_point(size = 1.5, alpha = 0.7) + xlim(-xlimit, xlimit) + geom_vline(xintercept = 2, colour = "blue", size = 1, alpha = 0.3) + geom_vline(xintercept = -2, colour = "blue", size = 1, alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) + xlab("Log 2 Fold-Change"))
	dev.off()
} #}}}
plotHeatmap <- function(y, heatmap_name, groups) { #  {{{
	    y <- y[,c(3,grep("sample", colnames(y)))]
	    y <- y[!duplicated(y[,1]),]
	    rownames(y) <- y[,1]
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
	    
	    png(paste(gsub("/",":", gsub(" ", "_", heatmap_name)), ".png", sep = ""), height = plot_height, width = plot_width, type = "cairo")
	    heatmap3(y[,-1], balanceColor = T, ColSideAnn=as.data.frame(groups), ColSideFun=function(x) showAnn(x), ColSideWidth = 0.8, main = heatmap_name, xlab = "Samples", ylab = "Genes", margins = c(15,10), cexRow = textsize)
	    dev.off()
} #}}}
plotDiffHeatmap <- function(y, groups, heatmap_name) { #  {{{
	    y <- y[,c(3,grep("sample", colnames(y)))]
	    y <- y[!duplicated(y[,1]),]
	    rownames(y) <- y[,1]
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

	    png(heatmap_name, height = plot_height, width = plot_width, type = "cairo")
	    heatmap3(y[,-1], balanceColor = T, ColSideAnn=as.data.frame(groups), ColSideFun=function(x) showAnn(x), ColSideWidth = 0.8, main = heatmap_name, xlab = "Samples", ylab = "Genes", margins = c(15,10), cexRow = textsize)
	    dev.off()
} #}}}
pathway <- function(x, species, keggDB, out.suffix) { #  {{{
	prev_dir <- getwd()
	out.suffix <- gsub(" ", "", out.suffix)
	new_dir <- paste(out.suffix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)
	exp.fc <- x$logFC
	names(exp.fc) <- x$entrezgene
	fc.kegg.p <- gage(exp.fc, gsets = keggDB, ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.kegg.p$greater)
	greater$Term_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.kegg.p$less)
	less$Term_ID <- trim(rownames(less))
	path_out <- merge(greater, less, by = "Term_ID")

	names(path_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	path_out$p.value <- ifelse(path_out$stat.mean.greater > 0, path_out$p.val.greater, path_out$p.val.less)
	path_out$FDR <- ifelse(path_out$stat.mean.greater > 0, path_out$q.val.greater, path_out$q.val.less)
	path_out <- path_out[,c(1,6,3,14,15)]
	names(path_out) <- c("Term_ID", "set.size", "mean_logFC", "P.value", "FDR")
	path_out <- path_out[!is.na(path_out$P.value),]
	row.names(path_out) <- path_out$Term_ID
	path_out <- path_out[order(path_out$P.value),]
	path_out$Genes <- "NULL"
	
	for (i in path_out$Term_ID) {
	path_entrezgenes <- kegg_db$kg.sets[[i]]
	path_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% path_entrezgenes, 3], collapse=","))
	}
	
	write.table(path_out, file = paste(out.suffix, "singleDirection_results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot barplots#
	try(plotSigTerms_KEGG_barplot(results = path_out, order = "top", num = 25, name = "Top25_Pvalue_significant_singleDirection_KEGG_barplot"))
	try(plotSigTerms_KEGG_barplot(results = path_out, order = "upregulated", num = 25, name = "Top25_upregulated_Pvalue_significant_singleDirection_KEGG_barplot"))
	try(plotSigTerms_KEGG_barplot(results = path_out, order = "downregulated", num = 25, name = "Top25_downregulated_Pvalue_significant_singleDirection_KEGG_barplot"))

	sel.greater <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "p.val"])
	path.ids.greater <- rownames(fc.kegg.p$greater)[sel.greater]
	sel.less <- fc.kegg.p$less[, "p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
	path.ids.less <- rownames(fc.kegg.p$less)[sel.less]

	fc.kegg.any <- gage(exp.fc, gsets = keggDB, ref = NULL, samp = NULL, same.dir = FALSE)
	path_out_any <- as.data.frame(fc.kegg.any$greater)
	path_out_any$Term_ID <- trim(rownames(path_out_any))
	path_out_any <- path_out_any[,c(7,5,2,3,4)]
	names(path_out_any) <- c("Term_ID", "set.size", "mean_logFC", "P.value", "FDR")
	path_out_any <- path_out_any[!is.na(path_out_any$P.value),]
	row.names(path_out_any) <- path_out_any$Term_ID
	path_out_any <- path_out_any[order(path_out_any$P.value),]
	path_out_any$Genes <- "NULL"
	
	for (i in path_out_any$Term_ID) {
	path_entrezgenes <- kegg_db$kg.sets[[i]]
	path_out_any[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% path_entrezgenes, 3], collapse=","))
	}

	sel.any <- fc.kegg.any$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.any$greater[, "p.val"])
	path.ids.any <- rownames(fc.kegg.any$greater)[sel.any]
	
	write.table(path_out_any, file = paste(out.suffix, "anyChange_results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot barplots#
	try(plotSigTerms_KEGG_barplot(results = path_out_any, order = "top", num = 25, name = "Top25_Pvalue_significant_anyChange_KEGG_barplot"))
	try(plotSigTerms_KEGG_barplot(results = path_out_any, order = "upregulated", num = 25, name = "Top25_upregulated_Pvalue_significant_anyChange_KEGG_barplot"))
	try(plotSigTerms_KEGG_barplot(results = path_out_any, order = "downregulated", num = 25, name = "Top25_downregulated_Pvalue_significant_anyChange_KEGG_barplot"))

	path.ids.all <- c(path.ids.greater, path.ids.less, path.ids.any)
	
	for (pid in path.ids.all) {
	    pathway.id <- substr(pid, 1, 8)
	    term_name <- gsub("/", ":", gsub(" ", "_", substr(pid,10, nchar(pid))))
	    try(pathview(gene.data =  exp.fc, pathway.id = pathway.id, species = species, out.suffix= paste(term_name, out.suffix, sep = "."), limit = list(gene=2, cpd=2), low = list(gene = "blue", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "orange", cpd = "orange")))
	    pid_entrezgenes <- kegg_db$kg.sets[[pid]]
	    pid_table <- x[x$entrezgene %in% pid_entrezgenes, ]
	    write.table(pid_table, file =  paste(pathway.id, term_name, out.suffix, "Differential_Expression_subset.xls", sep = "."), sep = "\t", quote = F, row.names = F)
	    file.remove(paste(pathway.id, ".png", sep = ""))
	    file.remove(paste(pathway.id, ".xml", sep = ""))
	     }
	setwd(prev_dir)
} #}}}
GO_BP <- function(x, counts, species, goDB, out.prefix = "GO_Biological_Process") { #  {{{
        prev_dir <- getwd()
	out.prefix <- gsub(" ", "", out.prefix)
	new_dir <- paste(out.prefix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)

	#GO Analysis#
	exp.fc <- x$logFC
	if ((biomart == "plants_mart") & (opt$feature == "gene")) {
	   names(exp.fc) <- x$Feature_ID
	}
	else if ((biomart == "plants_mart") & (opt$feature == "transcript")) {
	   names(exp.fc) <- x$ensembl_gene_id
	}else {
	   names(exp.fc) <- x$entrezgene
	}
	fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$BP], ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.go.p$greater)
	greater$Term_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.go.p$less)
	less$Term_ID <- trim(rownames(less))
	GO_out <- merge(greater, less, by = "Term_ID")
	names(GO_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, GO_out$p.val.less)
	GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, GO_out$q.val.less)
	GO_out <- GO_out[,c(1,6,3,14,15)]
	names(GO_out) <- c("Term_ID", "set.size", "mean_logFC", "P.value", "FDR")
	GO_out <- GO_out[!is.na(GO_out$P.value),]
	row.names(GO_out) <- GO_out$Term_ID
	GO_out <- GO_out[order(GO_out$P.value),]
	GO_out$Genes <- "NULL"
	
	for (i in GO_out$Term_ID) {
	GO_entrezgenes <- go_db$go.sets[[i]]
	GO_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% GO_entrezgenes, 3], collapse=","))
	GO_out
	}

	#Write Table
	write.table(GO_out, file = paste(out.prefix, "results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)
	
	#Plot barplots#
	try(plotSigTerms_GO_barplot(results = GO_out, order = "top", num = 25, name = "Top25_FDR_significant_GO_biological_processes_barplot"))
	try(plotSigTerms_GO_barplot(results = GO_out, order = "upregulated", num = 25, name = "Top25_upregulated_FDR_significant_GO_biological_processes_barplot"))
	try(plotSigTerms_GO_barplot(results = GO_out, order = "downregulated", num = 25, name = "Top25_downregulated_FDR_significant_GO_biological_processes_barplot"))

	#Plot Heatmaps#
	sel <- GO_out[GO_out$FDR <= 0.05, 1]
	for (gs in sel) {
	    #go.id <- substr(gs, 1, 10)
	    #outname = paste(paste(out.prefix, gsub("GO|:|/", "", gsub(" ","_", gs)), sep = "_"), ".png", sep = "")
	    #outname = paste(paste(out.prefix, gsub("GO|:|/", "", go.id), sep = "_"), ".png", sep = "")
	    
	    if ((biomart == "plants_mart") & (opt$feature == "gene")) {
	        z <- x[x$Feature_ID %in% go_db$go.sets[[gs]], ]
	    }
	    else if ((biomart == "plants_mart") & (opt$feature == "transcript")) {
	   	z <- x[x$Feature_ID %in% go_db$go.sets[[gs]], ]
	    } else {
	     	z <- x[x$entrezgene %in% go_db$go.sets[[gs]], ]
	    }

	    write.table(z, file = paste(gsub("/",":", gsub(" ", "_", gs)), "_Differential_Expression_subset.xls", sep = ""), sep = "\t", quote = F, row.names = F)
	    
	    if (nrow(z) > 3) {
	       hplot <- try(plotHeatmap(z, heatmap_name = gs, groups = groups))
	    }
	}
	setwd(prev_dir)
} #}}}
GO_MF <- function(x, counts, species, goDB, out.prefix = "GO_Molecular_Function") { #  {{{
    prev_dir <- getwd()
	out.prefix <- gsub(" ", "", out.prefix)
	new_dir <- paste(out.prefix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)

	#GO Analysis#
	exp.fc <- x$logFC
	if ((biomart == "plants_mart") & (opt$feature == "gene")) {
	   names(exp.fc) <- x$Feature_ID
	}
	else if ((biomart == "plants_mart") & (opt$feature == "transcript")) {
	   names(exp.fc) <- x$ensembl_gene_id
	}else {
	   names(exp.fc) <- x$entrezgene
	}
	fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$MF], ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.go.p$greater)
	greater$Term_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.go.p$less)
	less$Term_ID <- trim(rownames(less))
	GO_out <- merge(greater, less, by = "Term_ID")
	names(GO_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, GO_out$p.val.less)
	GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, GO_out$q.val.less)
	GO_out <- GO_out[,c(1,6,3,14,15)]
	names(GO_out) <- c("Term_ID", "set.size", "mean_logFC", "P.value", "FDR")
	GO_out <- GO_out[!is.na(GO_out$P.value),]
	row.names(GO_out) <- GO_out$Term_ID
	GO_out <- GO_out[order(GO_out$P.value),]
	GO_out$Genes <- "NULL"
	
	for (i in GO_out$Term_ID) {
	GO_entrezgenes <- go_db$go.sets[[i]]
	GO_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% GO_entrezgenes, 3], collapse=","))
	GO_out
	}

	#Write table#
	write.table(GO_out, file = paste(out.prefix, "results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot barplots#
	try(plotSigTerms_GO_barplot(results = GO_out, order = "top", num = 25, name = "Top25_FDR_significant_GO_molecular_functions_barplot"))
	try(plotSigTerms_GO_barplot(results = GO_out, order = "upregulated", num = 25, name = "Top25_upregulated_FDR_significant_GO_molecular_functions_barplot"))
	try(plotSigTerms_GO_barplot(results = GO_out, order = "downregulated", num = 25, name = "Top25_downregulated_FDR_significant_GO_molecular_functions_barplot"))
	
	#Plot Heatmaps#
	sel <- GO_out[GO_out$FDR <= 0.05, 1]
	for (gs in sel) {
	    #go.id <- substr(gs, 1, 10)
	    #outname = paste(paste(out.prefix, gsub("GO|:|/", "", gsub(" ","_", gs)), sep = "_"), ".png", sep = "")
	    #outname = paste(paste(out.prefix, gsub("GO|:|/", "", go.id), sep = "_"), ".png", sep = "")
	    
	    if ((biomart == "plants_mart") & (opt$feature == "gene")) {
	        z <- x[x$Feature_ID %in% go_db$go.sets[[gs]], ]
	    }
	    else if ((biomart == "plants_mart") & (opt$feature == "transcript")) {
	   	z <- x[x$ensembl_gene_id %in% go_db$go.sets[[gs]], ]
	    } else {
	     	z <- x[x$entrezgene %in% go_db$go.sets[[gs]], ]
	    }

	    write.table(z, file = paste(gsub("/",":", gsub(" ", "_", gs)), "_Differential_Expression_subset.xls", sep = ""), sep = "\t", quote = F, row.names = F)
	    
	    if (nrow(z) > 3){
	       hplot <- try(plotHeatmap(z, heatmap_name = gs, groups = groups))
	    }
	    
	}
	setwd(prev_dir)
} #}}}
plotGOHeatmap <- function(DIFF, GO_results,GO_number, groups) { #  {{{
	genelist <- strsplit(GO_results[grep(GO_number, GO_results$Term_ID),"Genes"], ",")
	heatmap_name <- GO_results[grep(GO_number, GO_results$Term_ID), "Term_ID"][[1]]
	y <- DIFF[DIFF$external_gene_name %in% genelist[[1]],]
	y <- y[,c(3,grep("sample", colnames(y)))]
        y <- y[!duplicated(y[,1]),]
        rownames(y) <- y[,1]
        plot_height <- min(c(max(c(1500, (nrow(y)*10))), 30000))
        plot_width <- min(c(max(c(1500, (nrow(y)*15))), 40000))
	textsize <- ifelse(nrow(y) <= 30, 1.8,1.4)
	tiff(paste(gsub("/",":", gsub(" ", "_", heatmap_name)), ".tiff", sep = ""), height = plot_height, width = plot_width)
  	heatmap3(y[,-1], balanceColor = T, ColSideAnn=as.data.frame(groups), ColSideFun=function(x) showAnn(x), ColSideWidth = 0.8, main = heatmap_name, margins = c(15,10), cexRow = textsize, cexCol = 2)
	dev.off()
}#}}}
plotGo_lineplot <-  function(DIFF, GO_results,GO_number, name) { #  {{{
	genelist <- strsplit(GO_results[grep(GO_number, GO_results$Term_ID),"Genes"], ",")
	lineplot_name <- GO_results[grep(GO_number, GO_results$Term_ID), "Term_ID"][[1]]
	y <- DIFF[DIFF$external_gene_name %in% genelist[[1]],]
	y <- y[,c("external_gene_name", "logFC", "CI.L", "CI.R", "P.Value")]
        y <- y[!duplicated(y$external_gene_name),]
	y <- y[order(y$P.Value),]
	y$Rank <- c(1:nrow(y))
	#tiff(file = name, height = 9, width = 16, units = 'in', dpi = 300)
	ggplot(y, aes(x = Rank)) + geom_ribbon(aes(ymin=CI.L, ymax=CI.R), alpha = 0.5) + geom_line(aes(y=logFC), colour = "red") + geom_line(aes(y = -log10(P.Value)), colour = "blue", size = 1) + geom_abline(intercept = 2, slope = 0) + xlab ("Significance Rank") + ylab("Log 2 fold-change (red) and -log10(P.Value) (blue)") + labs(title = lineplot_name) + theme(axis.text.y = element_text(hjust = 1, size = 12, colour = "black"))
	#dev.off()
}#}}}
plotGoPValue_lineplot <-  function(DIFF, GO_results,GO_number, name) { #  {{{
	genelist <- strsplit(GO_results[grep(GO_number, GO_results$Term_ID),"Genes"], ",")
	lineplot_name <- GO_results[grep(GO_number, GO_results$Term_ID), "Term_ID"][[1]]
	y <- DIFF[DIFF$external_gene_name %in% genelist[[1]],]
	y <- y[,c("external_gene_name", "logFC", "CI.L", "CI.R", "P.Value")]
        y <- y[!duplicated(y$external_gene_name),]
	y <- y[order(y$P.Value),]
	y$Rank <- c(1:nrow(y))
	#tiff(file = name, height = 9, width = 16, units = 'in', dpi = 300)
	ggplot(y, aes(x = Rank, y = -log10(P.Value))) + geom_line(colour = "blue", size = 1) + xlab ("Significance Rank") + ylab("-log10(P.Value)") + labs(title = lineplot_name) + theme(axis.text.y = element_text(hjust = 1, size = 12, colour = "black"))+ guides(color=FALSE)
	#dev.off()
}#}}}
plotGoFC_lineplot <-  function(DIFF, GO_results,GO_number, name) { #  {{{
	genelist <- strsplit(GO_results[grep(GO_number, GO_results$Term_ID),"Genes"], ",")
	lineplot_name <- GO_results[grep(GO_number, GO_results$Term_ID), "Term_ID"][[1]]
	y <- DIFF[DIFF$external_gene_name %in% genelist[[1]],]
	y <- y[,c("external_gene_name", "logFC", "CI.L", "CI.R", "P.Value")]
        y <- y[!duplicated(y$external_gene_name),]
	y <- y[order(-y$logFC),]
	y$Rank <- c(1:nrow(y))
	#tiff(file = name, height = 9, width = 16, units = 'in', dpi = 300)
	ggplot(y, aes(x = Rank)) + geom_ribbon(aes(ymin=CI.L, ymax=CI.R), alpha = 0.5) + geom_line(aes(y=logFC), colour = "red") + xlab ("Significance Rank") + geom_abline(intercept = 2, slope = 0) + ylab("Log 2 Fold-Change") + labs(title = lineplot_name) + theme(axis.text.y = element_text(hjust = 1, size = 12, colour = "black"))+ guides(color=FALSE)
	#dev.off()
}#}}}
plotSigTerms_GO_barplot <- function(results, order, num, name){ #  {{{
	if (order == "top") {
	   results <- results[c(1:num),]
	   }
	if (order == "upregulated") {
	   results <- results[results$FDR <= 0.05 & results$mean_logFC > 0,]
	   results <- results[order(-results$mean_logFC),]
	   results <- results[c(1:num),]
	   }
	if (order == "downregulated") {
	   results <- results[results$FDR <= 0.05 & results$mean_logFC < 0,]
	   results <- results[order(results$mean_logFC),]
	   results <- results[c(1:num),]
	   }
	   png(file = paste(name, "png", sep = "."), height = 1500, width = 2000)
	   print(ggplot(results, aes(y = -log10(P.value), x = Term_ID, fill = mean_logFC)) + geom_bar(stat = "identity") + coord_flip() + labs(title = name) +  scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + theme(axis.text.y = element_text(hjust = 1, size = 12, colour = "black")))
	   dev.off()
} #}}}
plotSigTerms_KEGG_barplot <- function(results, order, num, name){ #  {{{
	if (order == "top") {
	   results <- results[c(1:num),]
	   }
	if (order == "upregulated") {
	   results <- results[results$P.value <= 0.05 & results$mean_logFC > 0,]
	   results <- results[order(-results$mean_logFC),]
	   results <- results[c(1:num),]
	   }
	if (order == "downregulated") {
	   results <- results[results$P.value <= 0.05 & results$mean_logFC < 0,]
	   results <- results[order(results$mean_logFC),]
	   results <- results[c(1:num),]
	   }
	   png(file = paste(name, "png", sep = "."), height = 1500, width = 2000)
	   print(ggplot(results, aes(y = -log10(P.value), x = Term_ID, fill = mean_logFC)) + geom_bar(stat = "identity") + coord_flip() + labs(title = name) +  scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + theme(axis.text.y = element_text(hjust = 1, size = 12, colour = "black")))
	   dev.off()
}#}}}
enrichGOterms <- function(x, species = opt$species, out.prefix) {#  {{{
    require(clusterProfiler)
    if (species == "human") {
      orgdb <- "org.Hs.eg.db"
    } else if (species == "mouse") {
      orgdb <- "org.Mm.eg.db"
    } else if (species == "fly") {
      orgdb <- "org.Dm.eg.db"
    } else if (species == "rat") {
      orgdb <- "org.Rn.eg.db"
    } else if (species == "arabidopsis") {
      orgdb <- "org.At.tair.db"
    } else {
       	cat("Species not supported at this time, skipping enrichment analysis....\n")
    }

       	prev_dir <- getwd()
        out.prefix <- gsub(" ", "", out.prefix)
        new_dir <- paste(out.prefix, "GO_Enrichment_Results", sep = "_")
        dir.create(new_dir)
        setwd(new_dir)
    
	GO_BP_enrich <- enrichGO(gene = x$entrezgene, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "BP", readable = TRUE)
    	GO_BP_OUT <- as.data.frame(summary(GO_BP_enrich))
    	write.table(GO_BP_OUT, file = paste(out.prefix, "GO_Biological_Process_enrichment.xls", sep = "."), quote = F, sep = "\t", row.names = F)
    	png(file = paste(out.prefix, "top50_GO_Biological_Process_terms.png", sep = "."), height = 1500, width = 2000)
    	print(barplot(GO_BP_enrich, showCategory = 50)); dev.off()
    	png(file = paste(out.prefix, "GO_Biological_Process_Category_network_plot.png", sep = "."), height = 1500, width = 2000)
    	cnetplot(GO_BP_enrich, showCategory = 20, categorySize = "pvalue"); dev.off()

    	GO_MF_enrich <- enrichGO(gene = x$entrezgene, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "MF", readable = TRUE)
    	GO_MF_OUT <- as.data.frame(summary(GO_MF_enrich))
    	write.table(GO_MF_OUT, file = paste(out.prefix, "GO_Molecular_Function_enrichment.xls", sep = "."), quote = F, sep = "\t", row.names = F)
    	png(file = paste(out.prefix, "top50_GO_Molecular_Function_terms.png", sep = "."), height = 1500, width = 2000)
    	print(barplot(GO_MF_enrich, showCategory = 50)); dev.off()
    	png(file = paste(out.prefix, "GO_Molecular_Function_Category_network_plot.png", sep = "."), height = 1500, width = 2000)
    	cnetplot(GO_MF_enrich, showCategory = 20, categorySize = "pvalue"); dev.off()

	return(list(GO_BP_enrich, GO_MF_enrich))
	setwd(prev_dir)
 }#}}}
gseaGOterms <- function(x, species = opt$species, out.prefix) {#  {{{
    require(clusterProfiler)
    if (species == "human") {
      orgdb <- "org.Hs.eg.db"
    } else if (species == "mouse") {
      orgdb <- "org.Mm.eg.db"
    } else if (species == "fly") {
      orgdb <- "org.Dm.eg.db"
    } else if (species == "rat") {
      orgdb <- "org.Rn.eg.db"
    } else {
       	cat("Species not supported at this time, skipping enrichment analysis....\n")
    }

       	prev_dir <- getwd()
        out.prefix <- gsub(" ", "", out.prefix)
        new_dir <- paste(out.prefix, "GO_GSEA_Results", sep = "_")
        dir.create(new_dir)
        setwd(new_dir)

	gList <- x$logFC
	names(gList) <- x$entrezgene
	gList <- sort(gList, decreasing = TRUE)
    
	GO_BP_enrich <- gseGO(geneList = gList, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, ont = "BP")
    	GO_BP_OUT <- as.data.frame(summary(GO_BP_enrich))
    	write.table(GO_BP_OUT, file = paste(out.prefix, "GO_Biological_Process_GSEA_results.xls", sep = "."), quote = F, sep = "\t", row.names = F)

	if (nrow(GO_BP_OUT) >= 1) {
    	   #png(file = paste(out.prefix, "top50_GO_Biological_Process_terms.png", sep = "."), height = 1500, width = 2000)
    	   #print(barplot(GO_BP_enrich, showCategory = 50)); dev.off()
    	   #png(file = paste(out.prefix, "GO_Biological_Process_Category_network_plot.png", sep = "."), height = 1500, width = 2000)
    	   #cnetplot(GO_BP_enrich, showCategory = 20, categorySize = "pvalue", foldChange = gList); dev.off()
	   for (goterm in GO_BP_OUT$ID) {
	       goterm <- gsub("GO:", "", goterm)
	       png(file = paste(out.prefix, goterm, "GSEA_enrichment_plot.png", sep = "."), height = 1500, width = 2000)
	       gseaplot(GO_BP_enrich, geneSetID = goterm)
	       dev.off()
	   }
	}
	
	
    	GO_MF_enrich <- gseGO(geneList = gList, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, ont = "MF")
    	GO_MF_OUT <- as.data.frame(summary(GO_MF_enrich))
    	write.table(GO_MF_OUT, file = paste(out.prefix, "GO_Molecular_Function_GSEA_results.xls", sep = "."), quote = F, sep = "\t", row.names = F)
	
	if (nrow(GO_MF_OUT) >= 1) {
    	   #png(file = paste(out.prefix, "top50_GO_Molecular_Function_terms.png", sep = "."), height = 1500, width = 2000)
    	   #print(barplot(GO_MF_enrich, showCategory = 50)); dev.off()
    	   #png(file = paste(out.prefix, "GO_Molecular_Function_Category_network_plot.png", sep = "."), height = 1500, width = 2000)
    	   #cnetplot(GO_MF_enrich, showCategory = 20, categorySize = "pvalue", foldChange = gList); dev.off()
	   for (goterm in GO_MF_OUT$ID) {
	       goterm <- gsub("GO:", "", goterm)
	       png(file = paste(out.prefix, goterm, "GSEA_enrichment_plot.png", sep = "."), height = 1500, width = 2000)
	       gseaplot(GO_MF_enrich, geneSetID = goterm)
	       dev.off()
	   }
	}

    setwd(prev_dir)
}#}}}
gseaKEGGterms <- function(x, species = opt$species, out.prefix) { #  {{{
    require(clusterProfiler)
    if (species == "human") {
      orgdb <- "hsa"
    } else if (species == "mouse") {
      orgdb <- "mmu"
    } else if (species == "fly") {
      orgdb <- "dme"
    } else if (species == "rat") {
      orgdb <- "rno"
    } else {
       	cat("Species not supported at this time, skipping enrichment analysis....\n")
    }

       	prev_dir <- getwd()
        out.prefix <- gsub(" ", "", out.prefix)
        new_dir <- paste(out.prefix, "KEGG_GSEA_Results", sep = "_")
        dir.create(new_dir)
        setwd(new_dir)

	gList <- x$logFC
	names(gList) <- x$entrezgene
	gList <- sort(gList, decreasing = TRUE)
    
	KEGG_enrich <- gseKEGG(geneList = gList, organism = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05)
    	KEGG_OUT <- as.data.frame(summary(KEGG_enrich))
    	write.table(KEGG_OUT, file = paste(out.prefix, "KEGG_GSEA_results.xls", sep = "."), quote = F, sep = "\t", row.names = F)

	if (nrow(KEGG_OUT) >= 1) {
    	   #png(file = paste(out.prefix, "top50_KEGG_terms.png", sep = "."), height = 1500, width = 2000)
    	   #print(barplot(KEGG_enrich, showCategory = 50)); dev.off()
    	   #png(file = paste(out.prefix, "KEGG_Biological_Process_Category_network_plot.png", sep = "."), height = 1500, width = 2000)
    	   #cnetplot(KEGG_enrich, showCategory = 20, categorySize = "pvalue", foldChange = gList); dev.off()
	   for (KEGGterm in KEGG_OUT$ID) {
	       png(file = paste(out.prefix, KEGGterm, "GSEA_enrichment_plot.png", sep = "."), height = 1500, width = 2000)
	       gseaplot(KEGG_enrich, geneSetID = KEGGterm)
	       dev.off()
	   }
	}

   setwd(prev_dir)
}#}}}
merge_DIFFresults <- function(annotations, expression_data = NULL, write = TRUE) {#  {{{
	files <- list.files(pattern = "Differential_Expression.xls", recursive = T)
	x <- lapply(files, 'read.delim', header = T, stringsAsFactors = F)
	names(x) <- gsub("group_|_Differential_Expression.xls", "", basename(as.character(files)))
	avg_expression <- x[[1]][, c("Feature_ID", "AveExpr")]
	avg_expression <- avg_expression[!duplicated(avg_expression$Feature_ID),]
	for (Contrast in names(x)) {
    	x[[Contrast]] <- x[[Contrast]][,c("Feature_ID", "logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val")]
	    x[[Contrast]] <- x[[Contrast]][!duplicated(x[[Contrast]]$Feature_ID),]
	    x[[Contrast]]$Contrast <- c(as.character(Contrast))
	}

	require(reshape)
	merged_df <- merge_recurse(x)
	melted_df <- melt(merged_df)
	
	cast_df <- cast(melted_df, Feature_ID~Contrast + variable)
	cast_df<- merge(cast_df, avg_expression, by = "Feature_ID")
	if (!is.null(annotations)) {
		cast_df <- merge(annotations, cast_df, by = "Feature_ID")
	}
	if (!is.null(expression_data)) {
		OUT <- merge(cast_df, expression_data, by = "Feature_ID")
	} else {
		OUT <- cast_df
	}
	if (write == "TRUE") {
	   write.table(OUT, file = "Merged_differential_expression_results.xls", sep = "\t", quote = F, row.names = F)
	   melted_df
	} else {
	  OUT
	}
}#}}}
GLMcoefficient <- function(fit, coefficients, name, annotations, covariance = TRUE) { #  {{{
	 results_name <- gsub(" ", "", as.character(name))
	 dir.create(results_name)
         setwd(results_name)
	 
	 if (platform == "rnaseq") {
	   fit2 <- contrasts.fit(fit, coefficients=coefficients)
	   fit2 <- eBayes(fit2)
	 } else {
	   fit2 <- eBayes(fit2, trend = T, robust = T)
	 }

	 DIFF = topTable(fit2,n=Inf, sort = "none", confint = T)
	 if (!is.null(DIFF$logFC)) {
	    DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, (-1/(2^DIFF$logFC)))
	    
	    #Make an MA plot, showing genes that meet the FDR threshold in red
	    png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value)) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()
	 
	     plot_Volcano(DIFF, name = results_name)
	 }
	 
	  #Annotate results
	  DIFF$Feature_ID <- rownames(DIFF)
	  DIFF <- merge(DIFF, expression_data, by = "Feature_ID")
	  
	  if (!is.null(annotations)) {
	  	  DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		  DIFF <- DIFF[!duplicated(DIFF),]
	  }
	  
	    
 	  #Filter for significant results
	  if (covariance == "FALSE") {
	     DIFF$P.Value <- 1 - DIFF$P.Value
	     DIFF$adj.P.Val <- 1 - DIFF$adj.P.Val
	  }

	  DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
	  DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

	  #Write Data to a file
	    OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
	    OUT_DATA = format(OUT_DATA, round=5)
	    write.table(OUT_DATA, file=paste(results_name,"Differential_Expression.xls", sep = "_"),row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_SIG <- DIFF_SIG[order(DIFF_SIG$P.Value), ]
	    OUT_SIG = format(OUT_SIG, round=5)
	    write.table(OUT_SIG, file=paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_FDR <- DIFF_FDR[order(DIFF_FDR$adj.P.Val), ]
	    OUT_FDR = format(OUT_FDR, round=5)
	    write.table(OUT_FDR, file=paste(results_name, "FDR_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    if (nrow(DIFF_SIG) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	       	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	    }

	    if (nrow(DIFF_FDR) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	   }
	    if (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	   }
	   
	    if ((!is.null(kegg_db)) & (!is.null(DIFF$logFC))) {
	       pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "."))
	       if (go_species == "human" | go_species == "mouse") {
	       	  pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "."))
	       }
	   }       

	    if ((!is.null(go_db)) & (!is.null(DIFF$logFC))) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "."))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "."))
	    }

	    if ((!is.null(go_db)) & (is.null(DIFF$logFC))) {
	       enrichGOterms(DIFF_SIG, out.prefix = results_name)
	    }
        save.image(file = "DEX.results.RData")
	    
setwd(parent_dir)
}#}}}
GLMmatrix <- function(fit, matrix, annotations) {#  {{{
	 fit2 <- contrasts.fit(fit, matrix)
	 if (platform == "rnaseq") {
	   fit2 <- eBayes(fit2)
	 } else {
	   fit2 <- eBayes(fit2, trend = T, robust = T)
	 }

	for (test in colnames(fit2$coefficients)) {
	    results_name <- gsub(" ", "", as.character(test))
	    dir.create(results_name)
            setwd(results_name)
	    DIFF = topTable(fit2,sort = "none", n=Inf, confint = T, coef = test)
	    DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, (-1/(2^DIFF$logFC)))
	    #DIFF <- DIFF[,c(9,1:8)]

	    #Make an MA plot, showing genes that meet the FDR threshold in red
	    png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value)) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()
		
	    plot_Volcano(DIFF, name = results_name)
	    
	     #Annotate results
	     DIFF$Feature_ID <- rownames(DIFF)
	     DIFF <- merge(DIFF, expression_data, by = "Feature_ID")
	     if (!is.null(annotations)) {
	     	DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		DIFF <- DIFF[!duplicated(DIFF),]	  
 	    }

 
	    #Filter for significant results
	    DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
	    DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

	    #Write Data to a file
	    OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
	    OUT_DATA = format(OUT_DATA, round = 5)
	    write.table(OUT_DATA, file=paste(results_name,"Differential_Expression.xls", sep = "_"),row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_SIG <- DIFF_SIG[order(DIFF_SIG$P.Value), ]
	    OUT_SIG = format(OUT_SIG, round=5)
	    write.table(OUT_SIG, file=paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_FDR <- DIFF_FDR[order(DIFF_FDR$adj.P.Val), ]
	    OUT_FDR = format(OUT_FDR, round = 5)
	    write.table(OUT_FDR, file=paste(results_name, "FDR_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    if (nrow(DIFF_SIG) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	       	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	    }

	    if (nrow(DIFF_FDR) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	   }
	    if (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	   }
	   
	    if (!is.null(kegg_db)) {
	       	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "."))
	    	if (go_species == "human" | go_species == "mouse") {
	 	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "."))
	   	}	
	    }

	    if (!is.null(go_db)) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "."))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "."))
	    }

	    #enrichGOterms(DIFF_SIG, out.prefix = results_name)
        save.image(file = "DEX.results.RData")
	    
	    setwd(parent_dir)
	}

	#Create Venn Diagrams
	if (ncol(fit2$coefficients) <= 5) {
	results <- decideTests(fit2, adjust.method = "none")
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results);dev.off()
	results2 <- decideTests(fit2, adjust.method = "none", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results2);dev.off()
	results3 <- decideTests(fit2, adjust.method = "BH")
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results3);dev.off()
	results4 <- decideTests(fit2, adjust.method = "BH", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results4);dev.off()
	}

	#melted_results <- merge_DIFFresults(annotations = annotations, expression_data = expression_data, write = TRUE)
}#}}}
GLMmatrix <- function(fit, matrix, annotations) { # linearFC {{{
	 fit2 <- contrasts.fit(fit, matrix)
	 if (platform == "rnaseq") {
	   fit2 <- eBayes(fit2)
	 } else {
	   fit2 <- eBayes(fit2, trend = T, robust = T)
	 }

	for (test in colnames(fit2$coefficients)) {
	    results_name <- gsub(" ", "", as.character(test))
	    dir.create(results_name)
            setwd(results_name)
	    DIFF = topTable(fit2,sort = "none", n=Inf, confint = T, coef = test)
	    #DIFF <- DIFF[,c(9,1:8)]

	    #Make an MA plot, showing genes that meet the FDR threshold in red
	    png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value)) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()
		
	    plot_Volcano(DIFF, name = results_name)
	    
	     #Annotate results
	     DIFF$Feature_ID <- rownames(DIFF)
	    DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, (-1/(2^DIFF$logFC)))
        save(DIFF, fit2, expression_data, file = "DEX.results.RData")

	     DIFF <- merge(DIFF, expression_data, by = "Feature_ID")
	     if (!is.null(annotations)) {
	     	DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		DIFF <- DIFF[!duplicated(DIFF),]	  
 	    }

 
	    #Filter for significant results
	    DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
	    DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

	    #Write Data to a file
	    OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
	    OUT_DATA = format(OUT_DATA, round = 5)
	    write.table(OUT_DATA, file=paste(results_name,"Differential_Expression.xls", sep = "_"),row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_SIG <- DIFF_SIG[order(DIFF_SIG$P.Value), ]
	    OUT_SIG = format(OUT_SIG, round=5)
	    write.table(OUT_SIG, file=paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_FDR <- DIFF_FDR[order(DIFF_FDR$adj.P.Val), ]
	    OUT_FDR = format(OUT_FDR, round = 5)
	    write.table(OUT_FDR, file=paste(results_name, "FDR_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    if (nrow(DIFF_SIG) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	       	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	    }

	    if (nrow(DIFF_FDR) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	   }
	    if (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	   }
	   
	    if (!is.null(kegg_db)) {
	       	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "."))
	    	if (go_species == "human" | go_species == "mouse") {
	 	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "."))
	   	}	
	    }

	    if (!is.null(go_db)) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "."))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "."))
	    }

	    #enrichGOterms(DIFF_SIG, out.prefix = results_name)
	    
	    setwd(parent_dir)
	}

	#Create Venn Diagrams
	if (ncol(fit2$coefficients) <= 5) {
	results <- decideTests(fit2, adjust.method = "none")
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results);dev.off()
	results2 <- decideTests(fit2, adjust.method = "none", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results2);dev.off()
	results3 <- decideTests(fit2, adjust.method = "BH")
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results3);dev.off()
	results4 <- decideTests(fit2, adjust.method = "BH", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results4);dev.off()
	}

	#melted_results <- merge_DIFFresults(annotations = annotations, expression_data = expression_data, write = TRUE)
} #}}}
merge_DIFFresults <- function (annotations, expression_data = NULL, write = TRUE)# linearFC {{{
{
    files <- list.files(pattern = "Differential_Expression.xls",
                        recursive = T)
    x <- lapply(files, "read.delim", header = T, stringsAsFactors = F)
    names(x) <- gsub("group_|_Differential_Expression.xls", "",
                     basename(as.character(files)))
    avg_expression <- x[[1]][, c("Feature_ID", "AveExpr")]
    avg_expression <- avg_expression[!duplicated(avg_expression$Feature_ID),
                                     ]
    for (Contrast in names(x)) {
        x[[Contrast]] <- x[[Contrast]][, c("Feature_ID", "logFC", "linearFC",
                                           "CI.L", "CI.R", "P.Value", "adj.P.Val")]
                                                                               x[[Contrast]] <- x[[Contrast]][!duplicated(x[[Contrast]]$Feature_ID),
                                                                                                              ]
                                                                               x[[Contrast]]$Contrast <- c(as.character(Contrast))
    }
    require(reshape)
    merged_df <- merge_recurse(x)
    melted_df <- melt(merged_df)
    cast_df <- cast(melted_df, Feature_ID ~ Contrast + variable)
    cast_df <- merge(cast_df, avg_expression, by = "Feature_ID")
    if (!is.null(annotations)) {
        cast_df <- merge(annotations, cast_df, by = "Feature_ID")
    }
    if (!is.null(expression_data)) {
        OUT <- merge(cast_df, expression_data, by = "Feature_ID")
    }
    else {
        OUT <- cast_df
    }
    if (write == "TRUE") {
        write.table(OUT, file = "Merged_differential_expression_results.xls",
                    sep = "\t", quote = F, row.names = F)
        melted_df
    }
    else {
        OUT
    }
}#}}}
visAnt_module <- function(module, ngenes = "all", TOM = TOM, moduleColors = moduleColors, data = data) {#  {{{
        require(WGCNA)
       # Select module probes
       probes = names(data)
       inModule = (moduleColors==module);
       modProbes = probes[inModule];
       modAnnotations <- annotations[annotations$Feature_ID %in% modProbes, c("Feature_ID", "external_gene_name")]
       modAnnotations <- modAnnotations[!duplicated(modAnnotations$Feature_ID),]

       # Select the corresponding Topological Overlap
       modTOM = TOM[inModule, inModule];
       dimnames(modTOM) = list(modProbes, modProbes)
       modTOM <- modTOM[modAnnotations$Feature_ID, modAnnotations$Feature_ID]

       # Export the network into an edge list file VisANT can read
       if (ngenes == "all") {
              vis = exportNetworkToVisANT(modTOM, file = paste("VisANTInput-", module, ".txt", sep=""), weighted = TRUE, threshold = 0, probeToGene = modAnnotations)
	} else {
       	       nTop = ngenes;
       	       IMConn = softConnectivity(data[, modAnnotations$Feature_ID]);
       	       top = (rank(-IMConn) <= nTop)
       	       vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-", ngenes, "genes.txt", sep=""), weighted = TRUE, threshold = 0, probeToGene = modAnnotations)
	}
}#}}}
visAnt_FDR <- function(FDR, ngenes, TOM = TOM, data = data) {#  {{{
	   require(WGCNA)
	   FDR <- FDR$Feature_ID
	   # Select module probes
	   probes = names(data)
	   inFDR = (probes %in% FDR);
	   FDRprobes = probes[inFDR];
	   # Select the corresponding Topological Overlap
	   fdrTOM = TOM[inFDR, inFDR];
	   dimnames(fdrTOM) = list(FDRprobes, FDRprobes)

	   # Export the network into an edge list file VisANT can read
	    if (ngenes == "all") {
	       vis = exportNetworkToVisANT(fdrTOM, file = "VisANTInput-FDR.txt", weighted = TRUE, threshold = 0, probeToGene = annotations[!duplicated(annotations[,c("Feature_ID","external_gene_name")]),c("Feature_ID","external_gene_name")])
	     } else {
	       nTop = ngenes;
	       IMConn = softConnectivity(data[, FDRprobes]);
	       top = (rank(-IMConn) <= nTop)
	       vis = exportNetworkToVisANT(fdrTOM[top, top], file = paste("VisANTInput-", "FDR_topConnectivity", "-", ngenes, "genes.txt", sep=""), weighted = TRUE, threshold = 0, probeToGene = data.frame(annotations$Feature_ID, annotations$external_gene_name))
	       vis = exportNetworkToVisANT(fdrTOM[top, top],
	       file = "VisANTInput-FDR-top30.txt", weighted = TRUE, threshold = 0, probeToGene = data.frame(annotations$Feature_ID, annotations$external_gene_name))
	     }	     
} #}}}
moduleHeatmap <- function(module, Colors = moduleColors, data = data, groups) {#  {{{
       probes = names(data)
       inModule = (Colors==module);
       modProbes = probes[inModule];
       modExprs = t(data[,modProbes])
       require(heatmap3)

       plot_height <- min(c(max(c(1500, (nrow(modExprs)*10))), 20000))
       plot_width <- min(c(max(c(1500, (nrow(modExprs)*15))), 20000))

       if (nrow(modExprs) <= 50) {
	        textsize <- 1
	} else if ((nrow(modExprs) > 50) & (nrow(modExprs) <= 1500)) {
	    	textsize <- 0.8
	} else {
	      	textsize <- 0.1
	}

       modExprs <- as.data.frame(modExprs)
       modExprs$Feature_ID <- rownames(modExprs)
       modExprs <- merge(annotations[,c("Feature_ID", "external_gene_name")], modExprs, by = "Feature_ID")
       modExprs <- modExprs[,-1]
       modExprs <- modExprs[!duplicated(modExprs$external_gene_name),]
       rownames(modExprs) <- modExprs$external_gene_name
       
       tiff(file = paste(module, "Heatmap.tiff", sep = "_"), height = 1200, width = 2000)
       heatmap3(modExprs[,-1], balanceColor = T, ColSideAnn=as.data.frame(groups), ColSideFun=function(x) showAnn(x), ColSideWidth = 0.8, main = module, xlab = "Samples", ylab = "Genes", margins = c(15,10), cexRow = textsize)
       dev.off()
}

fdrHeatmap <-  function(FDR, data = data) {
       probes = names(data)
       FDR <- FDR$Feature_ID
       inFDR = (probes %in% FDR);
       FDRprobes = probes[inFDR];
       fdrExprs = t(data[,FDRprobes])
       require(heatmap3)
       tiff(file = "FDR_heatmap.tiff", height = 1200, width = 2000);heatmap3(fdrExprs, balanceColor = T, main = "FDR");dev.off()
}#}}}
moduleGO <- function(module, species = species, geneInfoCorr = geneInfoCorr, out.prefix) {#  {{{
     require(clusterProfiler)
     if (species == "human") {
      orgdb <- "org.Hs.eg.db"
     } else if (species == "mouse") {
      orgdb <- "org.Mm.eg.db"
     } else if (species == "fly") {
      orgdb <- "org.Dm.eg.db"
     } else if (species == "rat") {
      orgdb <- "org.Rn.eg.db"
     } else if (species == "worm") {
      orgdb <- "org.Ce.eg.db"
     } else if (species == "arabidopsis") {
      orgdb <- "org.At.tair.db"
     } else {
       	cat("Species not supported at this time, skipping enrichment analysis....\n")
     }
	EZ <- geneInfoCorr[geneInfoCorr$moduleColor == module & !is.na(geneInfoCorr$entrezgene),]
	EZ <- EZ[order(-EZ$Correlation.ModuleMembership),]
	GO_BP_enrich <- enrichGO(gene = EZ$entrezgene, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "BP", readable = TRUE)
    	GO_BP_OUT <- as.data.frame(GO_BP_enrich)
    	write.table(GO_BP_OUT, file = paste(out.prefix, "_GO_Biological_Process_enrichment.xls", sep = ""), quote = F, sep = "\t", row.names = F)
    	png(file = paste(out.prefix, "_top50_GO_Biological_Process_terms.png", sep = ""), height = 1500, width = 2000)
    	print(barplot(GO_BP_enrich, showCategory = 50))
	dev.off()
    	png(file = paste(out.prefix, "_GO_Biological_Process_Category_network_plot.png", sep = ""), height = 1500, width = 2000)
    	cnetplot(GO_BP_enrich, showCategory = 20, categorySize = "pvalue")
	dev.off()
	#png(file = paste(out.prefix, "_GO_Biological_Process_mapped_network_plot.png", sep = ""), height = 1500, width = 2000)
    	#emapplot(GO_BP_enrich, showCategory = 20, color = "pvalue"); dev.off()

    	GO_MF_enrich <- enrichGO(gene = EZ$entrezgene, OrgDb = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "MF", readable = TRUE)
    	GO_MF_OUT <- as.data.frame(GO_MF_enrich)
    	write.table(GO_MF_OUT, file = paste(out.prefix, "_GO_Molecular_Function_enrichment.xls", sep = ""), quote = F, sep = "\t", row.names = F)
    	png(file = paste(out.prefix, "_top50_GO_Molecular_Function_terms.png", sep = ""), height = 1500, width = 2000)
    	print(barplot(GO_MF_enrich, showCategory = 50))
	dev.off()
    	png(file = paste(out.prefix, "_GO_Molecular_Function_Category_network_plot.png", sep = ""), height = 1500, width = 2000)
    	cnetplot(GO_MF_enrich, showCategory = 20, categorySize = "pvalue")
	dev.off()
	#png(file = paste(out.prefix, "_GO_Molecular_Function_mapped_network_plot.png", sep = ""), height = 1500, width = 2000)
    	#emapplot(GO_MF_enrich, showCategory = 20, color = "pvalue"); dev.off()
}#}}}
moduleKEGG <- function(module, species = species, geneInfoCorr = geneInfoCorr, out.prefix) {#  {{{
     require(clusterProfiler)
     if (species == "human") {
      orgdb <- "hsa"
     } else if (species == "mouse") {
      orgdb <- "mmu"
     } else if (species == "fly") {
      orgdb <- "dme"
     } else if (species == "rat") {
      orgdb <- "rno"
     } else if (species == "worm") {
      orgdb <- "cel"
     } else if (species == "arabidopsis") {
      orgdb <- "ath"
     } else {
       	cat("Species not supported at this time, skipping enrichment analysis....\n")
     }

     if (species == "arabidopsis") {
     	EZ <- geneInfoCorr[geneInfoCorr$moduleColor == module & !is.na(geneInfoCorr$Feature_ID),]
	EZ <- EZ[order(-EZ$Correlation.ModuleMembership),]
	KEGG_enrich <- enrichKEGG(gene = EZ$Feature_ID, organism = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
     } else {
	EZ <- geneInfoCorr[geneInfoCorr$moduleColor == module & !is.na(geneInfoCorr$entrezgene),]
	EZ <- EZ[order(-EZ$Correlation.ModuleMembership),]
	KEGG_enrich <- enrichKEGG(gene = EZ$entrezgene, organism = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
     }
     
     KEGG_OUT <- as.data.frame(KEGG_enrich)
     write.table(KEGG_OUT, file = paste(out.prefix, "_KEGG_enrichment.xls", sep = ""), quote = F, sep = "\t", row.names = F)
     png(file = paste(out.prefix, "_top50_KEGG_terms.png", sep = ""), height = 1500, width = 2000)
     print(barplot(KEGG_enrich, showCategory = 50))
     dev.off()
     png(file = paste(out.prefix, "_KEGG_Category_network_plot.png", sep = ""), height = 1500, width = 2000)
     cnetplot(KEGG_enrich, showCategory = 20, categorySize = "pvalue")
     dev.off()
     
     #png(file = paste(out.prefix, "_KEGG_mapped_network_plot.png", sep = ""), height = 1500, width = 2000)
     #emapplot(KEGG_enrich, showCategory = 20, color = "pvalue"); dev.off()
     #KEGG_M_enrich <- enrichMKEGG(gene = EZ$entrezgene, organism = orgdb, maxGSSize = 1500, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
     #KEGG_M_OUT <- as.data.frame(KEGG_M_enrich)
     #write.table(KEGG_M_OUT, file = paste(out.prefix, "_KEGG_Module_enrichment.xls", sep = ""), quote = F, sep = "\t", row.names = F)
     #png(file = paste(out.prefix, "_top50_KEGG_Module_terms.png", sep = ""), height = 1500, width = 2000)
     #print(barplot(KEGG_M_enrich, showCategory = 50)); dev.off()
     #png(file = paste(out.prefix, "_KEGG_Module_Category_network_plot.png", sep = ""), height = 1500, width = 2000)
     #cnetplot(KEGG_M_enrich, showCategory = 20, categorySize = "pvalue"); dev.off()
     #png(file = paste(out.prefix, "_KEGG_Module_mapped_network_plot.png", sep = ""), height = 1500, width = 2000)
     #emapplot(KEGG_M_enrich, showCategory = 20, color = "pvalue"); dev.off()
}
#}}}
runWGCNA <- function(wgcna_data = datExpr, traits = datTraits, numSamples = nSamples, gene_annotations = annotations, wgcna_species = opt$species, threshold, correlationType = "signedHybrid") {#  {{{
	 require(WGCNA)
	 require(reshape)
	 require(ggplot2)
	 
	 cat("Merging differential expression results...\n")
	 #merged_results <- merge_DIFFresults(annotations = annotations, write = FALSE)
	 if (file.exists("Merged_differential_expression_results.xls")) {
	    merged_results <- read.delim(file = "Merged_differential_expression_results.xls", sep = "\t", stringsAsFactors = F)
	 } else {
	    merged_results <- read.delim(file = list.files(pattern = "Differential_Expression.xls", recursive = T), sep = "\t", header = T, stringsAsFactors = F)
	 }

	 dir.create("WGCNA")
	 setwd("WGCNA")
   	 cat("Creating topology overlap matrix...\n")

	 if (correlationType == "signedHybrid") {
	    TomType <- "signed"
   	    NetworkType <- "signed hybrid"
	 }
	 if (correlationType == "unsigned") {
   	    TomType <- "unsigned"
   	    NetworkType <- "unsigned"
	 }
	 if (correlationType == "signed") {
   	    TomType <- "signed"
   	    NetworkType <- "signed"
	 }

   	 net <- blockwiseModules(wgcna_data, maxBlockSize = nGenes, power = threshold, TOMType = TomType, networkType = NetworkType, minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, loadTOM = TRUE, saveTOMs = TRUE, saveTOMFileBase = "TopologyOverlapMatrix", verbose = 3)

   	 moduleLabels = net$colors
   	 moduleColors = labels2colors(net$colors)

   	 #Write table of detected modules and numbers of genes per module
   	 write.table(table(moduleColors), file = "Detected_modules_and_Number_of_Genes_per_Module.txt", sep = "\t", quote = F, row.names = F)

   	 #Plot dendrogram and modules
   	 tiff(file = "Sample_Dendrogram_and_Detected_Modules.tiff", height = 1000, width = 1500);plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05);dev.off()

   	 MEs = net$MEs;
   	 geneTree = net$dendrograms[[1]]

   	 ##Module Eigengene Significance##
   	 cat("Calculating module eigengene significance...\n")
   	 MEs0 = moduleEigengenes(wgcna_data, moduleColors)$eigengenes
   	 MEs = orderMEs(MEs0)
   	 modNames = substring(names(MEs), 3)
   	 moduleTraitCor = cor(MEs, traits, use = "p")
   	 moduleTraitPvalue = corPvalueStudent(moduleTraitCor, numSamples)

   	 #Create dataframe of Module Eigengene Significance and Correlations
   	 moduleTraitCor_melt <- melt.matrix(moduleTraitCor)
   	 names(moduleTraitCor_melt) <- c("Module_Color", "Trait", "Correlation")
   	 moduleTraitPvalue_melt <- melt.matrix(moduleTraitPvalue)
   	 names(moduleTraitPvalue_melt) <- c("Module_Color", "Trait", "P.value")
   	 moduleSignificance <- merge(moduleTraitCor_melt, moduleTraitPvalue_melt, by = c("Module_Color", "Trait"))
   	 moduleSignificance$Module_Color <- as.factor(moduleSignificance$Module_Color)
   	 moduleSignificance$Trait <- as.factor(moduleSignificance$Trait)
   	 moduleSignificance$Correlation <- as.numeric(format(as.numeric(moduleSignificance$Correlation), digits = 2))
   	 moduleSignificance$P.value <- as.numeric(format(as.numeric(moduleSignificance$P.value), digits = 2))

   	 #Create plot of correlations of module eigengenes to traits/phenotypes/conditions
   	 cat("Creating module and trait significance plot...\n")
   	 Sig_Matrix_Plot <- ggplot(moduleSignificance, aes(x = Trait, y = Module_Color, fill = Correlation)) + geom_tile() + geom_text(aes(Trait, Module_Color, label = paste(moduleSignificance$Correlation, moduleSignificance$P.value, sep = "\n")), size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, colour = "black")) + theme(axis.text.y = element_text(size = 20, colour = "black")) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", guide = "colourbar", na.value = "grey50") + labs(title = "Module and Trait Significance\n(Correlation and P-value)") + theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 22, face = "bold"), legend.title=element_text(size=20) , legend.text=element_text(size=18))
	 tiff(file = "Module_and_Trait_Eigengene_Significance_matrix.tiff", height = 4000, width = 5000);print(Sig_Matrix_Plot);dev.off()
   	 #ggsave(file = "Module_and_Trait_Eigengene_Significance_matrix.tiff", plot = Sig_Matrix_Plot);
   
	#Create dataframes of correlations of geneModuleMembership and Module Membership p-values
   	geneModuleMembership = as.data.frame(cor(wgcna_data, MEs, use = "p"))
   	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), numSamples))
   	names(geneModuleMembership) = sub("ME", "", names(geneModuleMembership))
   	names(MMPvalue) = sub("ME", "", names(MMPvalue))

   	#Melt geneModuleMembership and MMPvalue to prepare data for merging with results
   	geneModuleMembership$Feature_ID <- rownames(geneModuleMembership)
   	geneModuleMembership_melt <- melt(geneModuleMembership)
   	names(geneModuleMembership_melt) <- c("Feature_ID", "moduleColor", "Correlation.ModuleMembership")
   	geneModuleMembership_melt <- geneModuleMembership_melt[order(geneModuleMembership_melt$moduleColor, -geneModuleMembership_melt$Correlation.ModuleMembership),]
   	MMPvalue$Feature_ID <- rownames(MMPvalue)
   	MMPvalue_melt <- melt(MMPvalue)
   	names(MMPvalue_melt) <- c("Feature_ID", "moduleColor", "P.value.ModuleMembership")
   	MMPvalue_melt <- MMPvalue_melt[order(MMPvalue_melt$moduleColor, MMPvalue_melt$P.value.ModuleMembership),]

   	#Create dataframes of correlations and p-values of gene significance for each trait
   	geneTraitSignificance = as.data.frame(cor(wgcna_data, traits, use = "p"))
   	GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), numSamples))
   	names(geneTraitSignificance) = paste("Correlation.Gene.", names(traits), sep="")
   	names(GSPvalue) = paste("P.value.Gene.", names(traits), sep="")
      
	geneInfoPvalue0 <- data.frame(Feature_ID = names(wgcna_data), moduleColor = moduleColors, GSPvalue)
      	geneInfoPvalue0 <- merge(geneInfoPvalue0, geneModuleMembership_melt, by = c("Feature_ID", "moduleColor"))
      	geneInfoPvalue0 <- merge(geneInfoPvalue0, MMPvalue_melt, by = c("Feature_ID", "moduleColor"))
      	geneInfoPvalue <- geneInfoPvalue0[order(geneInfoPvalue0$moduleColor),]

   	geneInfoCorr0 <- data.frame(Feature_ID = names(wgcna_data), moduleColor = moduleColors, geneTraitSignificance)
   	geneInfoCorr0 <- merge(geneInfoCorr0, geneModuleMembership_melt, by = c("Feature_ID", "moduleColor"))
   	geneInfoCorr0 <- merge(geneInfoCorr0, MMPvalue_melt, by = c("Feature_ID", "moduleColor"))
   	geneInfoCorr <- geneInfoCorr0[order(geneInfoCorr0$moduleColor),]

	geneInfoDIFF0 <- data.frame(Feature_ID = names(wgcna_data), moduleColor = moduleColors)
   	geneInfoDIFF0 <- merge(geneInfoDIFF0, geneModuleMembership_melt, by = c("Feature_ID", "moduleColor"))
   	geneInfoDIFF0 <- merge(geneInfoDIFF0, MMPvalue_melt, by = c("Feature_ID", "moduleColor"))
	geneInfoDIFF <- merge(geneInfoDIFF0, merged_results, by = "Feature_ID")
	write.table(geneInfoDIFF, file = "WGCNA_Annotated_Differential_Expression_Results.xls", sep = "\t", quote = F, row.names = F)
	
   
	#Prepare and merge gene_annotations if available
	if (!is.null(gene_annotations)) {
   	   names(gene_annotations)[1] <- c("Feature_ID")
   	   geneInfoPvalue <- merge(gene_annotations, geneInfoPvalue, by = "Feature_ID")
   	   geneInfoPvalue$entrezgene <- as.character(geneInfoPvalue$entrezgene)
   	   geneInfoPvalue <- geneInfoPvalue[order(geneInfoPvalue$moduleColor),]
   	   geneInfoCorr <- merge(gene_annotations, geneInfoCorr, by = "Feature_ID")
   	   geneInfoCorr$entrezgene <- as.character(geneInfoCorr$entrezgene)
   	   geneInfoCorr <- geneInfoCorr[order(geneInfoCorr$moduleColor),]
   	   geneModuleMembershipCorr <- merge(gene_annotations, geneModuleMembership, by = "Feature_ID")
   	   geneModuleMembershipCorr$entrezgene <- as.character(geneModuleMembershipCorr$entrezgene)
   	   geneModuleMembershipPvalue <- merge(gene_annotations, MMPvalue, by = "Feature_ID")
   	   geneModuleMembershipPvalue$entrezgene <- as.character(geneModuleMembershipPvalue$entrezgene)
	} else {
	   geneModuleMembershipCorr <- geneModuleMembership[,c(ncol(geneModuleMembership),1:(ncol(geneModuleMembership) - 1))]
	   geneModuleMembershipPvalue <- MMPvalue[,c(ncol(MMPvalue), 1:(ncol(MMPvalue) - 1))]
	}
	
   	cat("Writing tables.....\n")
   	write.table(geneInfoPvalue, file = "WGCNA_module_and_trait_Pvalue_results.xls", sep = "\t", quote = F, row.names = F)
   	write.table(geneInfoCorr, file = "WGCNA_module_and_trait_Correlation_results.xls", sep = "\t", quote = F, row.names = F)
   	write.table(geneModuleMembershipCorr, file = "WGCNA_ModuleMembership_Correlation_results.xls", sep = "\t", quote = F, row.names = F)
   	write.table(geneModuleMembershipPvalue, file = "WGCNA_ModuleMembership_Pvalue_results.xls", sep = "\t", quote = F, row.names = F)


   	cat("Clustering, network identification, and significance testing complete, saving data.....\n\n")

	save(list = ls(all.names = TRUE), file = "WGCNA_analysis.RData") 


   	TOM = TOMsimilarityFromExpr(wgcna_data, power = threshold);

	
   	#Run GO analysis on modules if gene_annotations are available
	#if (!is.null(gene_annotations)) {
      	#    probes <- names(wgcna_data)
      	#    probe2annot <- match(probes,gene_annotations[,1])
      	#    allLLIDs = gene_annotations$entrezgene[probe2annot]
      	#    GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = wgcna_species, nBestP = 10);
      	#    GOtable = GOenr$bestPTerms[[4]]$enrichment
      	#    write.table(GOtable, file="GOenrichment_in_Modules.xls", sep = "\t", quote = F)
   	#}

	for (mod in modNames) {
	    prev.dir <- getwd()
	    dir.create(mod)
	    setwd(mod)
	    try(moduleHeatmap(module = mod, Colors = moduleColors, data = wgcna_data, groups = groups))
	    if (!is.null(gene_annotations)) {
	  	  try(visAnt_module(module = mod, ngenes = "all", TOM = TOM, moduleColors = moduleColors, data = wgcna_data))
		  graphics.off()
	  	  try(moduleGO(module = mod, species = wgcna_species, geneInfoCorr = geneInfoCorr, out.prefix = mod))
		  graphics.off()
		  try(moduleKEGG(module = mod, species = wgcna_species, geneInfoCorr = geneInfoCorr, out.prefix = mod))
		  graphics.off()
	    }
	    setwd(prev.dir)
	 }

   	cat("Analysis of module networks complete, saving data....\n\n")
	save(list = ls(all.names = TRUE), file = "WGCNA_analysis.RData")
	
	setwd(parent_dir)

	cat("FINISHED\n")
}#}}}

#####################################
##Work##
#####################################
##Read in and parse Pair File##
#####################################
#  {{{
cat("Reading in pair files....\n")
targets_file <- read.delim(file = samplekey, sep = "\t", stringsAsFactors = FALSE, header = F)
targets_file <- targets_file[,colSums(is.na(targets_file))<nrow(targets_file)]
parent_dir <- getwd()
#}}}

#####################################
##Statistical Models##
#####################################
if (opt$platform == "rnaseq") { #  {{{

	if (model == "standard") {
	   	if (opt$paired == "TRUE") {
		   	names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
		targets_file$group_ <- relevel(targets_file$group_, ref="control")
		design <- model.matrix(~group_, data = targets_file)
		null_design <- model.matrix(~1, data = targets_file)
	}
	if (model == "additive") {
	   	if (opt$paired == "TRUE") {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "factor_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "factor_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
		targets_file$group_ <- relevel(targets_file$group_, ref="control")
		targets_file$factor_ <- as.factor(targets_file$factor_)
		design <- model.matrix(~group_ + factor_, data = targets_file)
		null_design <- model.matrix(~factor_, data = targets_file)
	}
	if (model == "interaction") {
	   	if (opt$paired == "TRUE") {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "factor_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "factor_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
		targets_file$group_ <- relevel(targets_file$group_, ref = "control")
		targets_file$factor_ <- as.factor(targets_file$factor_)
		targets_file$factor_ <- relevel(targets_file$factor_, ref = "control")
		design <- model.matrix(~group_ + factor_ + group_:factor_, data = targets_file)
		null_design <- model.matrix(~1, data = targets_file)
	}
	if (model == "time-series") {
	   	if (opt$paired == "TRUE") {
		   	names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "time_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "time_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
	   	targets_file$group_ <- relevel(targets_file$group_, ref="control")
	   	targets_file$time_ <- as.factor(targets_file$time_)
	   	targets_file$time_ <- relevel(targets_file$time_, ref = "Day_0")
		#targets_file$subject_ <- as.factor(targets_file$subject_)
		#design <- model.matrix(~group_*time_, data = targets_file)
		design <- model.matrix(~group_ + time_ + group_:time_, data=targets_file)
		null_design <- model.matrix(~1, data=targets_file)
	}
	
	if (model == "fit") {
	
		if (opt$paired == "TRUE") {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
		design <- model.matrix(~0 + group_, data = targets_file)
		null_design <- model.matrix(~1, data = targets_file)
	}

	if (model == "fit-additive") {
	
		if (opt$paired == "TRUE") {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "batch_", "subject_")
			targets_file$subject_ <- as.factor(targets_file$subject_)
		} else {
			names(targets_file) <- c("index_name", "index", "sample", "path", "group_", "batch_")
		}
		targets_file$group_ <- fixnames(targets_file$group_)
		targets_file$group_ <- as.factor(targets_file$group_)
		targets_file$batch_ <- as.factor(targets_file$batch_)
		design <- model.matrix(~0 + group_ + batch_, data = targets_file)
		null_design <- model.matrix(~batch_, data = targets_file)
	}

	if (opt$feature == "gene") {
		count_file = "gene_counts.txt"
	} else {
		count_file = "sailfish/quant_bias_corrected.sf"
	}
	files <- ifelse(targets_file$path == "-", file.path(dirname(getwd()), targets_file$index, count_file), file.path("/scratch/gtac/analysis/rna_seq", targets_file$path, targets_file$index, count_file))
    #}}}
} else {#   {{{

	if (model == "standard") {
	   names(targets_file) <- c("files", "sample", "path", "group_")
	   targets_file$group_ <- fixnames(targets_file$group_)
	   targets_file$group_ <- as.factor(targets_file$group_)
	   targets_file$group_ <- relevel(targets_file$group_, ref="control")
	   design <- model.matrix(~group_, data = targets_file)
	   null_design <- model.matrix(~1, data = targets_file)
	}
	if (model == "additive") {
	   names(targets_file) <- c("files", "sample", "path", "group_", "sex_", "bmi_", "age_")
	   targets_file$group_ <- fixnames(targets_file$group_)
	   targets_file$group_ <- as.factor(targets_file$group_)
	   targets_file$sex_ <- as.factor(targets_file$sex_)
	   design <- model.matrix(~group_ + sex_ + bmi_ + age_, data = targets_file)
	   null_design <- model.matrix(~sex_ + bmi_ + age_, data = targets_file)
	}

	if (model == "interaction") {
	   names(targets_file) <- c("files", "sample", "path", "group_", "factor_")
	   targets_file$group_ <- fixnames(targets_file$group_)
	   targets_file$group_ <- as.factor(targets_file$group_)
	   targets_file$group_ <- relevel(targets_file$group_, ref="control")
	   targets_file$factor_ <- as.factor(targets_file$factor_)
	   targets_file$factor_ <- relevel(targets_file$factor_, ref="control")
	   design <- model.matrix(~group_*factor_, data = targets_file)
	   null_design <- model.matrix(~1, data = targets_file)
	}

	if (model == "time-series") {
	   names(targets_file) <- c("files", "sample", "path", "group_", "time_")
	   targets_file$group_ <- fixnames(targets_file$group_)
	   targets_file$group_ <- as.factor(targets_file$group_)
	   targets_file$group_ <- relevel(targets_file$group_, ref="control")
	   targets_file$time_ <- as.factor(targets_file$time_)
	   targets_file$time_ <- relevel(targets_file$time_, ref = "0hr")
	   design <- model.matrix(~group_*time_, data = targets_file)
	   #design <- model.matrix(~group_ + group_:time_, data=targets_file)
	   null_design <- model.matrix(~factor_, data = targets_file)
	}

	if (model == "fit") {
	   names(targets_file) <- c("files", "sample", "path", "group_", "sex_", "bmi_", "age_")
	   targets_file$group_ <- fixnames(targets_file$group_)
	   targets_file$group_ <- as.factor(targets_file$group_)
	   targets_file$sex_ <- as.factor(targets_file$sex_)
	   design <- model.matrix(~0 + group_ + sex_ + bmi_ + age_, data = targets_file)
	   null_design <- model.matrix(~sex_ + bmi_ + age_, data = targets_file)

	   #names(targets_file) <- c("files", "sample", "path", "group_")
	   #targets_file$group_ <- fixnames(targets_file$group_)
	   #targets_file$group_ <- as.factor(targets_file$group_)
	   #design <- model.matrix(~0+group_, data = targets_file)
	   #null_design <- model.matrix(~1_, data = targets_file)
	}
	
	files <- targets_file$files <- ifelse(targets_file$path == "-", file.path(dirname(getwd()), targets_file$files), file.path(targets_file$path, targets_file$files))
}#}}}

###################################
##Read in Data and BioMart Annotations##
#####################################
#  {{{
cat("Reading in raw data for normalization, model fitting, and downloading annotations...\n")
samples <- paste("sample", fixnames(targets_file$sample), sep = ".")
samples <- make.unique(samples)
groups <- fixnames(targets_file$group_)
#}}}

if (opt$feature == "gene") {#  {{{
   biomart_attributes <- c("ensembl_gene_id", "entrezgene", "external_gene_name", "gene_biotype", "external_gene_source", "transcript_count", "description")
   ensembl_filter <- "ensembl_gene_id"
}#}}}

if ((biomart == "plants_mart") & (opt$feature == "transcript")) {#  {{{
    biomart_attributes <- c("ensembl_transcript_id", "entrezgene", "ensembl_transcript_id", "transcript_biotype", "ensembl_gene_id", "external_gene_name", "description")
    ensembl_filter <- "ensembl_transcript_id"
}#}}}

if ((biomart != "plants_mart") & (opt$feature == "transcript")) {#  {{{
    biomart_attributes <- c("ensembl_transcript_id", "entrezgene", "external_transcript_name", "transcript_biotype", "ensembl_gene_id", "external_gene_name", "description")
    ensembl_filter <- "ensembl_transcript_id"
}#}}}

if (!is.null(reference_database)) {#  {{{
   mart = useMart(biomart, dataset = reference_database, host=host)
   #mart = useMart(biomart, dataset = reference_database)
}#}}}

if (platform == "agilent") {#  {{{
   data <- read.maimages(files = files, names = samples, source = "generic", columns = c("R" = "rMedianSignal", "Rb" = "rBGMedianSignal"), annotation = c("FeatureNum", "ProbeUID","ControlType","ProbeName", "GeneName", "SystematicName", "Description"))
   
   if (chip_type == "NULL") {
	   annotations <- read.delim(file = opt$custom_gene_annotations, sep = "\t", stringsAsFactors = F, header = T)
  } else if (chip_type == "efg_agilent_sureprint_g3_ge_8x60k_v3") {
    	   require("GEOquery")
    	   a <- getGEO("GPL21061")
 	   annotations <- Table(a)
	   names(annotations)[1] <- "Feature_ID"
  } else {
    	   annotations <- getBM(attributes = c(chip_type, biomart_attributes), filters = chip_type, values = data$genes$ProbeName, mart = mart)
  }

   data <- backgroundCorrect(data,method="normexp")
   d0 <- normalizeBetweenArrays(data,method="quantile")
   neg95 <- apply(d0$E[d0$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))  #Determine 95% of neg ctrl probes per array
   cutoff <- matrix(1.1*neg95,nrow(d0),ncol(d0),byrow=TRUE) 	#Filter for probes with intensities greater than 10% above neg95 in at least min
   isexpr <- rowSums(d0$E > cutoff) >= min(table(groups))
   d <- d0[d0$genes$ControlType==0 & isexpr, ]
   colnames(d$genes)[3] <- "Feature_ID"
   d <- avereps(d,ID=d$genes[,"Feature_ID"])
   d <- d[d$genes$Feature_ID %in% annotations[,1],] #Only keep annotated genes

   #Perform SVA
   if (opt$estLatentFactors == "TRUE") {
      sva.effects = sva(d$E,design,null_design)
      colnames(sva.effects$sv) <- paste("latentFactor_", 1:sva.effects$n.sv, "_", sep = "")
      design <- cbind(design, sva.effects$sv)
   }
}
#}}}
if (platform == "illumina") {#  {{{
   data <- read.idat(files, bgx)
   d <- neqc(data)
   rownames(d) <- d$genes$Probe_Id
   colnames(d$E) <- samples
   annotations <- getBM(attributes = c(chip_type, biomart_attributes), filters = chip_type, values = d$genes$Probe_Id, mart = mart)

   #Perform SVA
   if (opt$estLatentFactors == "TRUE") {
      sva.effects = sva(d$E,design,null_design)
      colnames(sva.effects$sv) <- paste("latentFactor_", 1:sva.effects$n.sv, "_", sep = "")
      design <- cbind(design, sva.effects$sv)
   }
}
#}}}
if (platform == "affy") {#  {{{
   data <- ReadAffy(filenames = files, sampleNames = samples)
   d <- rma(data)
   if (chip_type == "NULL") {
	   	annotations <- read.delim(file = opt$custom_gene_annotations, sep = "\t", stringsAsFactors = F, header = T)
   	} else {
   		annotations <- getBM(attributes = c(chip_type, biomart_attributes), filters = chip_type, values = rownames(exprs(d)), mart = mart)
	}
   #Perform SVA
   if (opt$estLatentFactors == "TRUE") {
      sva.effects = sva(exprs(d),design,null_design)
      colnames(sva.effects$sv) <- paste("latentFactor_", 1:sva.effects$n.sv, "_", sep = "")
      design <- cbind(design, sva.effects$sv)
   }
}
#}}}
if (platform == "rnaseq") {#  {{{
   targets <- data.frame("files" = files, "group" = groups, "description" = samples)
   if (opt$feature == "gene") {
      data <- readDGE(targets, header = T, columns = c(1,7), comment.char = "#")
   } else {
      data <- readDGE(targets, header = F, columns = c(1,7), comment.char = "#")
   }
   colnames(data) <- c(samples)

   if (!is.null(reference_database)) {
   	annotations <- getBM(attributes = biomart_attributes, filters = ensembl_filter, values = rownames(data), mart = mart)
   } else if ((!is.null(opt$custom_gene_annotations)) & (opt$feature == "gene")) {
     	annotations <- read.delim(file = opt$custom_gene_annotations, sep = "\t", stringsAsFactors = F)
   } else if ((!is.null(opt$custom_transcript_annotations)) & (opt$feature == "transcript")) {
     	annotations <- read.delim(file = opt$custom_transcript_annotations, sep = "\t", stringsAsFactors = F)
   } else {
   annotations <- data.frame("Feature_ID" = rownames(data))
   }

   #Filter out ribosomal and lowly expressed genes
   cat("Filtering lowly expressed genes and annotated ribosomal genes....\n")
   if (!is.null(annotations)) {
      	ribosomal_genes <- annotations[grep(pattern = "rRNA", annotations[,grep(pattern = "biotype", names(annotations))]),1]
   	data <- data[!rownames(data) %in% ribosomal_genes,]
   }

   if (opt$miRNA == "TRUE") {
      miRNA_genes <- annotations[grep(pattern = "miRNA", annotations[,grep(pattern = "biotype", names(annotations))]),1]
      data <- data[rownames(data) %in% miRNA_genes,]
   }

   if (opt$transgenic == "TRUE") {
      transgenic_genes <- grep(pattern = "transgenic", rownames(data))
      data <- data[transgenic_genes,]
   }

   cpm.d <- cpm(data)
   
   #CPM_THRESHOLD  = 1
   #MIN_REP_NUMBER  = min(table(groups))
   #EXPRESSED = (rowSums(cpm.d > CPM_THRESHOLD) >= MIN_REP_NUMBER)

    library(readxl)
    keep.genes <-  read_excel("/scratch/gtac/analysis/rna_seq/8998_3_magee-s4483/code/July\ 2019\ RNAseq\ FLT3ITD\ and\ ME-Lin28\ expressed\ genes.xlsx", 2)[[1]]
    #keep.chipid <- annotations[match(keep.genes, annotations$external_gene_name), "ensembl_gene_id"]
    keep.chipid <- keep.genes
   #EXPRESSED = (rowSums(cpm.d > 0) >= 1) & (rownames(data) %in% keep.chipid)
   EXPRESSED = (rownames(data) %in% keep.chipid)
   
   data <- data[EXPRESSED, ]
   
   cat("Final number of features and samples of filtered dataset:\n")
   print(dim(data))
   data$samples$lib.size <- colSums(data$counts)
   cat("Final library size for samples of filtered_dataset:\n")
   print(data$samples$lib.size)
   cat("\n\n")
   
   #Rescale counts by library size and VOOM transform counts
   data <- calcNormFactors(data, method = "TMM")
   
   #Perform SVA
   if (opt$estLatentFactors == "TRUE") {
      sva.effects = svaseq(data$counts,design,null_design)
      if (!is.null(dim(sva.effects$sv))) {
      	 colnames(sva.effects$sv) <- paste("latentFactor_", 1:sva.effects$n.sv, "_", sep = "")
	 design <- cbind(design, sva.effects$sv)
      } else {
      latentFactor <- sva.effects$sv
      design <- cbind(design, latentFactor)
      }
   }
   
   #Perform Voom transformation and Mean-Variance Modeling
   cat("Fitting for feature level mean-variance weights and sample quality weights with Voom.....\n")
   pdf(file = "Voom_Variance_Model_Fit_and_Sample_Weights.pdf", paper = "a4", title = "Voom Mean-Variance Model Fit and Sample Weights")
   d <- voomWithQualityWeights(data, design, plot = TRUE)
 
   if (opt$paired == "TRUE") {
       cat("Fitting duplicate correlation for paired samples...\n")
       corfit <- duplicateCorrelation(d,design,block=targets_file$subject_)
       d <- voomWithQualityWeights(data, design, plot=TRUE, block=targets_file$subject_, correlation=corfit$consensus)
   }
   dev.off()
   
   chip_type <- ensembl_filter
   d$genes <- data.frame("Feature_ID" = rownames(d))
}
#}}}
names(annotations)[1] <- c("Feature_ID")

 if (platform == "affy") {#  {{{
    expression_data <- as.data.frame(exprs(d))#}}}
 } else {#  {{{
    expression_data <- as.data.frame(d$E)
}#}}}
expression_data$Feature_ID <- rownames(expression_data)

cat("Data read in and filtered, normalized for statistical analysis, and annotations downloaded, proceeding with QC....\n")
#####################################
##Make an MDS plot##
#####################################
#  {{{
cat("Making MDS plot...\n")
png(file="MDS.png", height=600, width=600)						
MDS = plotMDS(d, main="MDS Plot", labels=samples) 
ggMDS_coord <- data.frame(x = MDS$x, y = MDS$y, condition = groups)
print(ggplot(ggMDS_coord, aes(x = x, y = y, colour = condition, label = rownames(ggMDS_coord), position = "dodge")) + geom_point(size = 6) + geom_text(hjust = 0, vjust = 2, size = 4) + geom_point(size = 6) + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("Multi-Dimensional Scaling Plot") + xlim(min(MDS$x)-3,max(MDS$x)+3) + ylim(min(MDS$y)-3,max(MDS$y)+3))
dev.off()

if (opt$estLatentFactors == "TRUE" & opt$model != "additive") {#  {{{
	if (opt$paired == "TRUE") {
		d2 <- removeBatchEffect(d, batch = targets_file$subject_, covariates = sva.effects$sv)
	} else {
		d2 <- removeBatchEffect(d, covariates = sva.effects$sv)	
	}
	png(file="MDS_batchEffectRemoved.png", height=600, width=600)				
	MDS_noBatch = plotMDS(d2, main="MDS Plot", labels=samples) 
	ggMDS_noBatch_coord <- data.frame(x = MDS_noBatch$x, y = MDS_noBatch$y, condition = groups)
	print(ggplot(ggMDS_noBatch_coord, aes(x = x, y = y, colour = condition, label = rownames(ggMDS_noBatch_coord), position = "dodge")) + geom_point(size = 6) + geom_text(hjust = 0, vjust = 2, size = 4) + geom_point(size = 6) + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("Multi-Dimensional Scaling Plot: Batch Effects Removed") + xlim(min(MDS$x)-3,max(MDS$x)+3) + ylim(min(MDS$y)-3,max(MDS$y)+3))
	dev.off()
}#}}}
#}}}

if (opt$paired == "TRUE" & opt$model != "additive" & opt$estLatentFactors == "FALSE") {#  {{{

	d2 <- removeBatchEffect(d, batch = targets_file$subject_)

	png(file="MDS_batchEffectRemoved.png", height=600, width=600)				
	MDS_noBatch = plotMDS(d2, main="MDS Plot", labels=samples) 
	ggMDS_noBatch_coord <- data.frame(x = MDS_noBatch$x, y = MDS_noBatch$y, condition = groups)
	print(ggplot(ggMDS_noBatch_coord, aes(x = x, y = y, colour = condition, label = rownames(ggMDS_noBatch_coord), position = "dodge")) + geom_point(size = 6) + geom_text(hjust = 0, vjust = 2, size = 4) + geom_point(size = 6) + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("Multi-Dimensional Scaling Plot: Batch Effects Removed") + xlim(min(MDS$x)-3,max(MDS$x)+3) + ylim(min(MDS$y)-3,max(MDS$y)+3))
	dev.off()
}#}}}

if (opt$model == "additive") {#  {{{
	if (opt$paired == "TRUE" & opt$estLatentFactors == "FALSE") {
		d2 <- removeBatchEffect(d, batch = targets_file$factor_, batch2 = targets_file$subject_)
	} else if (opt$paired == "TRUE" & opt$estLatentFactors == "TRUE") {
		d2 <- removeBatchEffect(d, batch = targets_file$factor_, batch2 = targets_file$subject_, covariates = sva.effects$sv)
	} else if (opt$paired == "FALSE" & opt$estLatentFactors == "TRUE") {
		d2 <- removeBatchEffect(d, batch = targets_file$factor_, covariates = sva.effects$sv)
	} else {
		d2 <- removeBatchEffect(d, batch = targets_file$factor_)
	}
	png(file="MDS_batchEffectRemoved.png", height=600, width=600)
	MDS_noBatch = plotMDS(d2, main="MDS Plot", labels=samples) 
	ggMDS_noBatch_coord <- data.frame(x = MDS_noBatch$x, y = MDS_noBatch$y, condition = groups)
	print(ggplot(ggMDS_noBatch_coord, aes(x = x, y = y, colour = condition, label = rownames(ggMDS_noBatch_coord), position = "dodge")) + geom_point(size = 6) + geom_text(hjust = 0, vjust = 2, size = 4) + geom_point(size = 6) + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("Multi-Dimensional Scaling Plot: Batch Effects Removed") + xlim(min(MDS$x)-3,max(MDS$x)+3) + ylim(min(MDS$y)-3,max(MDS$y)+3))
	dev.off()
}
#}}}
if (opt$model == "fit-additive") {#  {{{
	if (opt$paired == "TRUE" & opt$estLatentFactors == "FALSE") {
		d2 <- removeBatchEffect(d, batch = targets_file$batch_, batch2 = targets_file$subject_)
	} else if (opt$paired == "TRUE" & opt$estLatentFactors == "TRUE") {
		d2 <- removeBatchEffect(d, batch = targets_file$batch_, batch2 = targets_file$subject_, covariates = sva.effects$sv)
	} else if (opt$paired == "FALSE" & opt$estLatentFactors == "TRUE") {
		d2 <- removeBatchEffect(d, batch = targets_file$batch_, covariates = sva.effects$sv)
	} else {
		d2 <- removeBatchEffect(d, batch = targets_file$batch_)
	}
	png(file="MDS_batchEffectRemoved.png", height=600, width=600)
	MDS_noBatch = plotMDS(d2, main="MDS Plot", labels=samples) 
	ggMDS_noBatch_coord <- data.frame(x = MDS_noBatch$x, y = MDS_noBatch$y, condition = groups)
	print(ggplot(ggMDS_noBatch_coord, aes(x = x, y = y, colour = condition, label = rownames(ggMDS_noBatch_coord), position = "dodge")) + geom_point(size = 6) + geom_text(hjust = 0, vjust = 2, size = 4) + geom_point(size = 6) + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("Multi-Dimensional Scaling Plot: Batch Effects Removed") + xlim(min(MDS$x)-3,max(MDS$x)+3) + ylim(min(MDS$y)-3,max(MDS$y)+3))
	dev.off()
}
#}}}

#####################################
##Make a 3-D PCA plot##
#####################################
#  {{{
try(plot_PCA(expression_data, groups = groups))

if (exists("d2")) {
   try(plot_PCA(as.data.frame(d2), groups = groups, name = "PCA_plot_batchEffectRemoved"))
}#}}}

#####################################
##Make an Expression Box Plot##
#####################################
#  {{{
cat("Making a boxplot of filtered transformed expression values...\n")
expression_melt <- melt(expression_data)
names(expression_melt)[2] <- "Sample"
png(file = "Expression_boxplot.png", height = 1000, width = 1600)
ggplot(expression_melt, aes(x = Sample, y = value, colour = Sample)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black")) + theme(axis.text.y = element_text(size = 8)) + labs(title = "Box Plot of Voom Transformed Data")
dev.off()#}}}

#####################################
##Make a Correlation Plot##
#####################################
#  {{{
cat("Making Spearman correlation plot...\n")
if (platform == "affy") {
   correlations <- round(cor(exprs(d), use = "complete.obs", method = "spearman"), digits = 2)
} else {
   correlations <- round(cor(d$E, use = "complete.obs", method = "spearman"), digits = 2)
}
corr_melt <- melt(correlations)
names(corr_melt) <- c("Samples_1", "Samples_2", "Correlation")
png("Spearman_correlation_matrix.png", height = 1000, width = 1600)
ggplot(corr_melt, aes(x = Samples_1, y = Samples_2, fill = Correlation)) + geom_tile() + geom_text(aes(Samples_1, Samples_2, label = Correlation), size = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "grey50")) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradient(low="green", high="red") + labs(title = "Spearman Correlation Matrix")
dev.off()
#}}}

#####################################
##QC with WGCNA##
#####################################
#  {{{
cat("Loading WGCNA...\n")
require(DESeq2)
require(WGCNA)
allowWGCNAThreads()
if (opt$feature == "gene") {
	if (platform == "rnaseq") {
	   	datExpr <-  as.data.frame(t(varianceStabilizingTransformation(data$counts)))
	} else {
	    if (platform == "affy") {
	    datExpr <- as.data.frame(t(exprs(d)))
	    } else {
	    datExpr <- as.data.frame(t(d$E))
	    }	  
	}

	nGenes = ncol(datExpr)
	nSamples = nrow(datExpr)

	#Create datTraits
	datTraits <- as.data.frame(design)
	if (ncol(datTraits) == 2 & mean(datTraits[,1]) == 1) {
	      datTraits <- as.data.frame(model.matrix(~0+group_, data = targets_file))
	}
	if (ncol(datTraits) > 2 & mean(datTraits[,1]) == 1) {
	      datTraits <- datTraits[,-1]
	}

	##Cluster samples to detect outliers##
	sampleTree = hclust(dist(datExpr), method = "average")
	# Convert traits to a color representation: white means low, red means high, grey means missing entry
	traitColors = numbers2colors(datTraits, signed = FALSE);
	png(file = "Sample_and_Trait_Dendrogram.png", height = 1000, width = 1500);plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap", marAll = c(5,20,5,5), cex.colorLabels = 1, cex.dendroLabels = 1, cex.rowText = 1);dev.off()

   # Create topology and determine thresholds
   powers = c(c(1:10), seq(from = 12, to=30, by=2))
   sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
   sftThresholds <- data.frame("Threshold" = sft$fitIndices[,1], "Topology" = round(-sign(sft$fitIndices[,3])*sft$fitIndices[,2], digits = 3))
   
   if (max(sftThresholds$Topology) >= 0.9) {
      bestThreshold <- min(sftThresholds[sftThresholds$Topology > 0.9,1])
      cat("The best WGCNA power threshold parameter is:\n")
      cat(bestThreshold)
   } else if ((max(sftThresholds$Topology) <= 0.9) & (max(sftThresholds$Topology) >= 0.8)) {
      bestThreshold <- min(sftThresholds[sftThresholds$Topology > 0.8,1])
      cat("The best WGCNA power threshold parameter is:\n")
      cat(bestThreshold)
   } else {
      bestThreshold <- 9
      cat("There is no best WGCNA power threshold parameter\n")
      bestThreshold <- NULL
   }

   png(file = "Scale-Free_Topology_softThresholding_plot.png", height = 1000, width = 1500);par(mfrow = c(1,2));cex1 = 0.9;plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");abline(h=0.90,col="red");plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"));text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");dev.off()
}
#}}}

cat("QC complete, fitting statistical model....\n")
#####################################
##Fit and Test Statistical Models##
#####################################
##Fit Model###  {{{
if (model == "fit" | model == "fit-additive") {
   #Fit the model to the design of the experiment
   	if (platform == "rnaseq") {
	   if (opt$paired == "TRUE") {
   	      fit <- lmFit(d,design,block=targets_file$subject_,correlation=corfit$consensus)
	   } else {
  	     fit <- lmFit(d,design)
	   }	     
	} else {
	  fit <- lmFit(d, design)
	}

	if (opt$platform != "rnaseq") {
	   png(file="SA_mean-variance_plot.png", height = 600, width = 600);plotSA(fit, xlab="Average Log 2 Expression", ylab="Log 2 Residual Standard Deviation(sigma)", main="SA_mean-variance_Plot", cex = 0.5);dev.off()
	}
}#}}}
##Additive Model###  {{{
if (model == "additive") {
	
	#Fit the model to the design of the experiment
	if (platform == "rnaseq") {
	   if (opt$paired == "TRUE") {
   	      fit <- lmFit(d,design,block=targets_file$subject_,correlation=corfit$consensus)
   	      fit <- eBayes(fit)
	   } else {
  	     fit <- eBayes(lmFit(d,design))
	   }	     
	} else {
	  fit <- lmFit(d, design)
	  fit <- eBayes(fit, trend = T, robust = T)
	}

	if (opt$platform != "rnaseq") {
	   png(file="SA_mean-variance_plot.png", height = 600, width = 600);plotSA(fit, xlab="Average Log 2 Expression", ylab="Log 2 Residual Standard Deviation(sigma)", main="SA_mean-variance_Plot", cex = 0.5);dev.off()
	}

	DIFF = topTable(fit, sort = "none", n=Inf, confint = T, coef = 2)
	DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, (-1/(2^DIFF$logFC)))
	#DIFF <- DIFF[,c(9,10,1:8)]
	
	#Make an MA plot, showing genes that meet the FDR threshold in red
	png(file="MA_plot.png", height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value)) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()

	plot_Volcano(DIFF)
	
	#Annotate results
	DIFF$Feature_ID <- rownames(DIFF)
	DIFF <- merge(DIFF, expression_data, by = "Feature_ID")

	if (!is.null(annotations)) {
		DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		DIFF <- DIFF[!duplicated(DIFF),]
	}

	#Filter for significant results
	DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
	DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

	#Write Data to a file
	OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
	OUT_DATA = format(OUT_DATA, round=3)
	write.table(OUT_DATA, file="Differential_Expression.xls", row.names=FALSE, sep = "\t", quote = FALSE)

	OUT_SIG <- DIFF_SIG[order(DIFF_SIG$P.Value), ]
	OUT_SIG = format(OUT_SIG, round = 3)
	write.table(OUT_SIG, file="Unadjusted_Pvalue_Significant_Differential_Expression_Subset.xls", row.names=FALSE, sep = "\t", quote = FALSE)

	OUT_FDR <- DIFF_FDR[order(DIFF_FDR$adj.P.Val), ]
	OUT_FDR = format(OUT_FDR, round = 3)
	write.table(OUT_FDR, file="FDR_Significant_Differential_Expression_Subset.xls", row.names=FALSE, sep = "\t", quote = FALSE)

	    if (nrow(DIFF_SIG) > 3) {
	       plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png")
	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) > 3) {
	       plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png")
	    }

	    if (nrow(DIFF_FDR) > 3) {
	       plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = "FDR_Significant_Differential_Expression_Heatmap.png")
	   }
	    if  (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3) {
	       plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png")
	    }

	if (!is.null(kegg_db)) {
	   pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = "KEGG_Signaling_and_Metabolism")
	   if (go_species == "human" | go_species == "mouse") {
	 	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = "KEGG_Disease")
	   }	
	}		
	if (!is.null(go_db)) {
	   GO_BP(DIFF, species = go_species, goDB = go_db)
	   GO_MF(DIFF, species = go_species, goDB = go_db)
	}
	#enrichGOterms(DIFF_SIG, out.prefix = "Significant")

	##Run WGCNA##
	if (!is.null(bestThreshold)) {
		runWGCNA(threshold = bestThreshold)
	}
}#}}}
##Time-Series, Interaction, and Standard Models###  {{{
if (model == "time-series" | model == "interaction" | model == "standard") {

	#Fit the model to the design of the experiment
	if (platform == "rnaseq") {
	   if (opt$paired == "TRUE") {
   	      fit <- lmFit(d,design,block=targets_file$subject_,correlation=corfit$consensus)
   	      fit <- eBayes(fit)
	   } else {
  	     fit <- eBayes(lmFit(d,design))
	   }	     
	} else {
	  fit <- lmFit(d, design)
	  fit <- eBayes(fit, trend = T, robust = T)
	}
	
	if (opt$platform != "rnaseq") {
	   png(file="SA_mean-variance_plot.png", height = 600, width = 600);plotSA(fit, xlab="Average Log 2 Expression", ylab="Log 2 Residual Standard Deviation(sigma)", main="SA_mean-variance_Plot", cex = 0.5);dev.off()	
	}
	   
	for (test in colnames(fit$coefficients)[2:length(colnames(fit$coefficients))]) {
	    results_name <- gsub(" ", "", as.character(test))
	    dir.create(results_name)
            setwd(results_name)
	    DIFF = topTable(fit,sort = "none", n=Inf, confint = T, coef = test)
	    DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, (-1/(2^DIFF$logFC)))
	    #DIFF <- DIFF[,c(9,1:8)]	    

	    #Make an MA plot, showing genes that meet the FDR threshold in red
	    png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = P.Value)) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_gradient2(low = "red", mid = "blue", high = "grey50", midpoint = 0.4) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()
			
		plot_Volcano(DIFF, name = results_name)

	    
	    #Annotate results
	    DIFF$Feature_ID <- rownames(DIFF)
	    DIFF <- merge(DIFF, expression_data, by = "Feature_ID")

	    if (!is.null(annotations)) {
	     	DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		DIFF <- DIFF[!duplicated(DIFF),]
 	    }     
	    
	    #Filter for significant results
	    DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
	    DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

	    #Write Data to a file
	    OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
	    OUT_DATA = format(OUT_DATA, round = 5)
	    write.table(OUT_DATA, file=paste(results_name,"Differential_Expression.xls", sep = "_"),row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_SIG <- DIFF_SIG[order(DIFF_SIG$P.Value), ]
	    OUT_SIG = format(OUT_SIG, round=5)
	    write.table(OUT_SIG, file=paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    OUT_FDR <- DIFF_FDR[order(DIFF_FDR$adj.P.Val), ]
	    OUT_FDR = format(OUT_FDR, round = 5)
	    write.table(OUT_FDR, file=paste(results_name, "FDR_Significant_Differential_Expression_Subset.xls", sep = "_"), row.names=FALSE, sep = "\t", quote = FALSE)

	    if (nrow(DIFF_SIG) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	       	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	    }

	    if (nrow(DIFF_FDR) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_")))
	   }
	    if (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3) {
	       try(plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_")))
	   }
	       
	    if (!is.null(kegg_db)) {
	       	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "."))
	 	#pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "."))
	    }

	    if (!is.null(go_db)) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "."))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "."))
	    }

	    #enrichGOterms(DIFF_SIG, out.prefix = results_name)

	    setwd(parent_dir)
	}
	
	if (ncol(design) > 2) {
	   x <- merge_DIFFresults(annotations = annotations, expression_data = expression_data)
	}

	##Run WGCNA##
	if (!is.null(bestThreshold)) {
		runWGCNA(threshold = bestThreshold)
	}
}#}}}

#####################################
##Save Data##
#####################################
project_info <- sessionInfo()
save.image(file = "DEX_analysis.RData")

#####################################
##Report Errors##
#####################################
warnings()
cat("Done\n")
