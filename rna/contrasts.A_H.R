library(limma)
library(gage)
library(pathview)
library(ggplot2)
library(reshape)
library(heatmap3)

options(bitmapType = "cairo")

load("DEX_analysis.RData")

if(FALSE) {
pathway <- function(x, species, keggDB, out.suffix) { #  {{{
	prev_dir <- getwd()
	new_dir <- paste(out.suffix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)
	exp.fc <- x$logFC
	names(exp.fc) <- x$entrezgene
	fc.kegg.p <- gage(exp.fc, gsets = keggDB, ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.kegg.p$greater)
	greater$KEGG_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.kegg.p$less)
	less$KEGG_ID <- trim(rownames(less))
	path_out <- merge(greater, less, by = "KEGG_ID")

	names(path_out) <- c("KEGG_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	path_out$p.value <- ifelse(path_out$stat.mean.greater > 0, path_out$p.val.greater, path_out$p.val.less)
	path_out$FDR <- ifelse(path_out$stat.mean.greater > 0, path_out$q.val.greater, path_out$q.val.less)
	path_out <- path_out[,c(1,6,3,14,15)]
	names(path_out) <- c("KEGG_ID", "set.size", "mean_logFC", "P.value", "FDR")
	path_out <- path_out[!is.na(path_out$P.value),]
	row.names(path_out) <- path_out$KEGG_ID
	path_out <- path_out[order(path_out$P.value),]
	path_out$Genes <- "NULL"
	
	for (i in path_out$KEGG_ID) {
	path_entrezgenes <- kegg_db$kg.sets[[i]]
	path_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% path_entrezgenes, 3], collapse=","))
	path_out
	}
	
	write.table(path_out, file = paste(out.suffix, "singleDirection_results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)
	
	sel.greater <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "p.val"])
	path.ids.greater <- rownames(fc.kegg.p$greater)[sel.greater]
	sel.less <- fc.kegg.p$less[, "p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
	path.ids.less <- rownames(fc.kegg.p$less)[sel.less]

	fc.kegg.any <- gage(exp.fc, gsets = keggDB, ref = NULL, samp = NULL, same.dir = FALSE)
	path_out_any <- as.data.frame(fc.kegg.any$greater)
	path_out_any$KEGG_ID <- trim(rownames(path_out_any))
	path_out_any <- path_out_any[,c(7,5,2,3,4)]
	names(path_out_any) <- c("KEGG_ID", "set.size", "mean_logFC", "P.value", "FDR")
	path_out_any <- path_out_any[!is.na(path_out_any$P.value),]
	row.names(path_out_any) <- path_out_any$KEGG_ID
	path_out_any <- path_out_any[order(path_out_any$P.value),]
	path_out_any$Genes <- "NULL"
	
	for (i in path_out_any$KEGG_ID) {
	path_entrezgenes <- kegg_db$kg.sets[[i]]
	path_out_any[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% path_entrezgenes, 3], collapse=","))
	path_out_any
	}
	
	write.table(path_out_any, file = paste(out.suffix, "anyChange_results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)


	sel.any <- fc.kegg.any$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.any$greater[, "p.val"])
	path.ids.any <- rownames(fc.kegg.any$greater)[sel.any]

	path.ids.all <- substr(c(path.ids.greater, path.ids.less, path.ids.any), 1, 8)
	#pv.out.list <- sapply(path.ids.all, function(pid) pathview(gene.data =  exp.fc, pathway.id = pid, species = species, out.suffix=out.suffix, limit = list(gene=2, cpd=2), low = list(gene = "blue", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "orange", cpd = "orange")))
	setwd(prev_dir)
} #}}}
GO_CC <- function(x, counts, species, goDB, out.prefix = "GO_Cellular_Component") { #  {{{
        prev_dir <- getwd()
	new_dir <- paste(out.prefix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)

	#GO Analysis#
	exp.fc <- x$logFC
	names(exp.fc) <- x$entrezgene
	fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$CC], ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.go.p$greater)
	greater$GO_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.go.p$less)
	less$GO_ID <- trim(rownames(less))
	GO_out <- merge(greater, less, by = "GO_ID")
	names(GO_out) <- c("GO_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, GO_out$p.val.less)
	GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, GO_out$q.val.less)
	GO_out <- GO_out[,c(1,6,3,14,15)]
	names(GO_out) <- c("GO_ID", "set.size", "mean_logFC", "P.value", "FDR")
	GO_out <- GO_out[!is.na(GO_out$P.value),]
	row.names(GO_out) <- GO_out$GO_ID
	GO_out <- GO_out[order(GO_out$P.value),]
	GO_out$Genes <- "NULL"
	
	for (i in GO_out$GO_ID) {
	GO_entrezgenes <- go_db$go.sets[[i]]
	GO_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% GO_entrezgenes, 3], collapse=","))
	GO_out
	}
	
	write.table(GO_out, file = paste(out.prefix, "results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot Heatmaps#
	#sel <- GO_out[GO_out$P.value <= 0.05, 1]
	#for (gs in sel) {
		#go.id <- substr(gs, 1, 10)
		#outname = paste(paste(out.prefix, gsub("GO|:|/", "", go.id), sep = "_"), ".png", sep = "")
		#z <- x[x$entrezgene %in% go_db$go.sets[[gs]], ]
		
		#if ((nrow(z) <= 1500) & (nrow(z) > 3)){
		   #hplot <- try(plotHeatmap(z, heatmap_name = outname))
		#}
	#}
	setwd(prev_dir)
} #}}}
GO_BP <- function(x, counts, species, goDB, out.prefix = "GO_Biological_Process") { #  {{{
        prev_dir <- getwd()
	new_dir <- paste(out.prefix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)

	#GO Analysis#
	exp.fc <- x$logFC
	names(exp.fc) <- x$entrezgene
	fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$BP], ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.go.p$greater)
	greater$GO_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.go.p$less)
	less$GO_ID <- trim(rownames(less))
	GO_out <- merge(greater, less, by = "GO_ID")
	names(GO_out) <- c("GO_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, GO_out$p.val.less)
	GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, GO_out$q.val.less)
	GO_out <- GO_out[,c(1,6,3,14,15)]
	names(GO_out) <- c("GO_ID", "set.size", "mean_logFC", "P.value", "FDR")
	GO_out <- GO_out[!is.na(GO_out$P.value),]
	row.names(GO_out) <- GO_out$GO_ID
	GO_out <- GO_out[order(GO_out$P.value),]
	GO_out$Genes <- "NULL"
	
	for (i in GO_out$GO_ID) {
	GO_entrezgenes <- go_db$go.sets[[i]]
	GO_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% GO_entrezgenes, 3], collapse=","))
	GO_out
	}
	
	write.table(GO_out, file = paste(out.prefix, "results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot Heatmaps#
	#sel <- GO_out[GO_out$P.value <= 0.05, 1]
    #for (gs in sel) {
        #go.id <- substr(gs, 1, 10)
        #outname = paste(paste(out.prefix, gsub("GO|:|/", "", go.id), sep = "_"), ".png", sep = "")
        #z <- x[x$entrezgene %in% go_db$go.sets[[gs]], ]
        
        #if ((nrow(z) <= 1500) & (nrow(z) > 3)){
           #hplot <- try(plotHeatmap(z, heatmap_name = outname))
        #}
    #}
	setwd(prev_dir)
}#}}}
GO_MF <- function(x, counts, species, goDB, out.prefix = "GO_Molecular_Function") { #  {{{
        prev_dir <- getwd()
	new_dir <- paste(out.prefix, "Results", sep = "_")
	dir.create(new_dir)
        setwd(new_dir)

	#GO Analysis#
	exp.fc <- x$logFC
	names(exp.fc) <- x$entrezgene
	fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$MF], ref = NULL, samp = NULL)
	greater <- as.data.frame(fc.go.p$greater)
	greater$GO_ID <- trim(rownames(greater))
	less <- as.data.frame(fc.go.p$less)
	less$GO_ID <- trim(rownames(less))
	GO_out <- merge(greater, less, by = "GO_ID")
	names(GO_out) <- c("GO_ID", "p.geomean.greater", "stat.mean.greater", "p.val.greater", "q.val.greater", "set.size.greater", "exp1.greater", "p.geomean.less", "stat.mean.less", "p.val.less", "q.val.less", "set.size.less", "exp1.less")
	GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, GO_out$p.val.less)
	GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, GO_out$q.val.less)
	GO_out <- GO_out[,c(1,6,3,14,15)]
	names(GO_out) <- c("GO_ID", "set.size", "mean_logFC", "P.value", "FDR")
	GO_out <- GO_out[!is.na(GO_out$P.value),]
	row.names(GO_out) <- GO_out$GO_ID
	GO_out <- GO_out[order(GO_out$P.value),]
	GO_out$Genes <- "NULL"
	
	for (i in GO_out$GO_ID) {
	GO_entrezgenes <- go_db$go.sets[[i]]
	GO_out[i,"Genes"] <- as.character(paste(x[x$entrezgene %in% GO_entrezgenes, 3], collapse=","))
	GO_out
	}
	
	write.table(GO_out, file = paste(out.prefix, "results.xls", sep = "_"), sep = "\t", quote = F, row.names = F)

	#Plot Heatmaps#
	#sel <- GO_out[GO_out$P.value <= 0.05, 1]
    #for (gs in sel) {
        #go.id <- substr(gs, 1, 10)
        #outname = paste(paste(out.prefix, gsub("GO|:|/", "", go.id), sep = "_"), ".png", sep = "")
        #z <- x[x$entrezgene %in% go_db$go.sets[[gs]], ]

        #if ((nrow(z) <= 1500) & (nrow(z) > 3)){
           #hplot <- try(plotHeatmap(z, heatmap_name = outname))
        #}
        
    #}
	setwd(prev_dir)
} #}}}
GLMcontrast <- function(fit, test, name, annotations) { #  {{{
	 results_name <- as.character(name)
	 dir.create(results_name)
         setwd(results_name)

	 fit2 <- contrasts.fit(fit, coefficients=test)
	 #fit2 <- eBayes(fit2, trend = T, robust = T)
	 fit2 <- eBayes(fit2)
	 DIFF = topTable(fit2,n=Inf, sort = "none", confint = T)
	 DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, -1/(2^DIFF$logFC))
	 #DIFF <- DIFF[,c(9,10,1:8)]
	 DIFF <- DIFF[,c(1,10,2:9)]

	 #Make an MA plot, showing genes that meet the FDR threshold in red
	  png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = (adj.P.Val <= 0.05))) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()

	  #Annotate results
 	  if (platform == "affy") {
	   	   DIFF <- cbind(DIFF, as.data.frame(exprs(d)))
	   } else {
	   	   DIFF <- cbind(DIFF, as.data.frame(d$E))
	   }

	  DIFF$Feature_ID <- rownames(DIFF)

	  if (!is.null(annotations)) {
	  	  DIFF <- merge(annotations, DIFF, by = "Feature_ID")
		  DIFF <- DIFF[!duplicated(DIFF),]
	  }

 	  #Filter for significant results
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

	    if (nrow(DIFF_SIG) <= 500) {
	       plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_"))
	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) <= 1500) {
	       plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_"))
	    }

	    if ((nrow(DIFF_FDR) > 3) & (nrow(DIFF_FDR) <= 500)) {
	       plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_"))
	   }
	   if ((nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) <= 500) & (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3)) {
	      plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_"))
	   }

	    if (!is.null(kegg_db)) {
	       pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "_"))
	       if (go_species == "human" | go_species == "mouse") {
	       	       pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "_"))
	    }
	}       

	    if (!is.null(go_db)) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "_"))
	       GO_CC(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Cellular_Component", sep = "_"))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "_"))
	    }
setwd(parent_dir)
} #}}}
GLMmatrix <- function(fit, matrix, annotations) { #  {{{
	 fit2 <- contrasts.fit(fit, matrix)
	 #fit2 <- eBayes(fit2, trend = T, robust = T)
	 fit2 <- eBayes(fit2)

	for (test in colnames(fit2$coefficients)) {
	    results_name <- as.character(test)
	    dir.create(results_name)
            setwd(results_name)
	    DIFF = topTable(fit2,sort = "none", n=Inf, confint = T, coef = test)
	    DIFF$linearFC <- ifelse(DIFF$logFC > 0, 2^DIFF$logFC, -1/(2^DIFF$logFC))
	    DIFF <- DIFF[,c(9,10,1:8)]

	    #Make an MA plot, showing genes that meet the FDR threshold in red
	     png(file=paste(results_name, "MA_plot.png", sep = "_"), height=600, width=600);print(ggplot(DIFF, aes(x = AveExpr, y = logFC, colour = (adj.P.Val <= 0.05))) + geom_point(size = 4, alpha = 0.7) + geom_abline(intercept = -2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + geom_abline(intercept = 0, slope = 0, colour = "black", size = 1, alpha = 0.3) + geom_abline(intercept = 2, slope = 0, colour = "blue", size = 1, alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) + xlab("Average Log-Expression") + ylab("Log 2 Fold-change") + ggtitle("MA Plot"));dev.off()

	     #Annotate results
	        if (platform == "affy") {
	   	   DIFF <- cbind(DIFF, as.data.frame(exprs(d)))		   
		} else {
	   	   DIFF <- cbind(DIFF, as.data.frame(d$E))
		}

	     DIFF$Feature_ID <- rownames(DIFF)
	     
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

	    if (nrow(DIFF_SIG) <= 1500) {
	       plotDiffHeatmap(DIFF_SIG, groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap.png", sep = "_"))
	    }
	    if (nrow(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,]) <= 1500) {
	       plotDiffHeatmap(DIFF_SIG[DIFF_SIG$logFC >= 2 | DIFF_SIG$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "Unadjusted_Pvalue_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_"))
	    }

	    if ((nrow(DIFF_FDR) > 3) & (nrow(DIFF_FDR) <= 1500)) {
	       plotDiffHeatmap(DIFF_FDR, groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap.png", sep = "_"))
	   }
	   if ((nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) <= 1500) & (nrow(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,]) > 3)) {
	      plotDiffHeatmap(DIFF_FDR[DIFF_FDR$logFC >= 2 | DIFF_FDR$logFC <= -2,], groups = groups, heatmap_name = paste(results_name, "FDR_Significant_Differential_Expression_Heatmap_logFC_2.png", sep = "_"))
	   }
	   
	    if (!is.null(kegg_db)) {
	       	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$sigmet.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Signaling_and_Metabolism", sep = "_"))
	    	if (go_species == "human" | go_species == "mouse") {
	 	pathway(DIFF, keggDB = kegg_db$kg.set[kegg_db$dise.idx], species = kegg_species, out.suffix = paste(results_name, "KEGG_Disease", sep = "_"))
	   	}	
	    }

	    if (!is.null(go_db)) {
	       GO_BP(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Biological_Process", sep = "_"))
	       GO_CC(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Cellular_Component", sep = "_"))
	       GO_MF(DIFF, species = go_species, goDB = go_db, out.prefix = paste(results_name, "GO_Molecular_Function", sep = "_"))
	    }
	    setwd(parent_dir)
	}
	results <- decideTests(fit2, adjust.method = "none")
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results);dev.off()
	results2 <- decideTests(fit2, adjust.method = "none", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_p-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results2);dev.off()
	results3 <- decideTests(fit2, adjust.method = "BH")
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05.png", sep = "_"), height=1200, width=2000);vennDiagram(results3);dev.off()
	results4 <- decideTests(fit2, adjust.method = "BH", lfc = 2)
	png(file=paste(results_name, "Venn_Diagram_q-value_0.05_lfc_2.png", sep = "_"), height=1200, width=2000);vennDiagram(results4);dev.off()
} #}}}
#debug(GLMmatrix)
#plotDiffHeatmap <- function(...){}
#GO_BP <- GO_MF <- GO_CC <- pathway <- function(...){}
#contrast.matrix_1 <- makeContrasts("EvsF"=group_E - group_F,
                                   #levels = design)
}

GLMmatrix <- function(fit, matrix, annotations) { #  {{{
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
merge_DIFFresults <- function (annotations, expression_data = NULL, write = TRUE)#  {{{
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

contrast.matrix_1 <- makeContrasts("AvsC"=group_A - group_C,
                                   "AvsE"=group_A - group_E,
                                   "AvsG"=group_A - group_G,
                                   "CvsG"=group_C - group_G,
                                   "EvsG"=group_E - group_G,
                                   "BvsD"=group_B - group_D,
                                   "BvsF"=group_B - group_F,
                                   "BvsH"=group_B - group_H,
                                   "DvsH"=group_D - group_H,
                                   "FvsH"=group_F - group_H,
                                   levels = design)
#debug(GLMmatrix)
#plotDiffHeatmap <- function(...){}
#GO_BP <- GO_MF <- GO_CC <- pathway <- function(...){}
#contrast.matrix_1 <- makeContrasts("EvsF"=group_E - group_F,
                                   #levels = design)
GLMmatrix(fit, matrix = contrast.matrix_1, annotations)
merge_DIFFresults(expression_data = expression_data, annotations = annotations)
#runWGCNA(threshold = bestThreshold)
