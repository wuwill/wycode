library(Seurat)
library(readr)
library(stringr)
library(ggplot)
library(dplyr)

my.export.findmarkers <- function(x, path, ...) { #{{{
    library(readr)
    x <- data.frame(gene=rownames(x), x)
    write_tsv(x, path, ...)
} #}}}
my_QC <- function(dt, dt_name, antibodies = NULL, traits = NULL, suffix = ".beforeQC0", #  {{{
                  percent_pattern = list(mt = "^(MT|mt|Mt)-", #  {{{
                                         ribo = "^R(ps|pl|PS|PL)"),
                  plot_lim = list(nFeature_RNA = c(0, 2000),
                                  nCount_RNA = c(0, 5000),
                                  percent.mt = c(0, 2),
                                  percent.ribo = c(0, 2)),
                  qc_lim = list(nFeature_RNA = c(250, 4000),
                                  percent.mt = 5,
                                  percent.ribo = 5),
                  norm_method = "sct",
                  vars.to.regress = c("percent.mt", "percent.ribo"),
                  ...) { # ... used in data normalization }}}
    test.na <- function(x) length(x) == 1 && is.na(x)
    calc_percentage <- function(dt, patterns) {#{{{
        for(name in names(patterns))
            dt <- dt %>% PercentageFeatureSet(pattern = patterns[[name]],
                                              col.name = paste0("percent.", name))
        return(dt)
    }#}}}
    qc_plots0 <- function(dt, antibodies, lim = list(nFeature_RNA = c(0, 2000), #  {{{
                                                    nCount_RNA = c(0, 5000),
                                                    percent.mt = c(0, 2),
                                                    percent.ribo = c(0, 2)
                                                    )) {
        ploti <- function(feature) {#{{{
            print(VlnPlot(dt, features = feature, pt.size = 0.00) +
                  NoLegend() +
                  geom_boxplot(width=0.1,fill="white"))
            if(feature %in% names(lim)) {
                print(VlnPlot(dt, features = feature, pt.size = 0.01) + ylim(lim[[feature]]))
                print(VlnPlot(dt, features = feature, pt.size = 0) +
                      ylim(lim[[feature]]) +
                      geom_boxplot(width=0.1,fill="white"))
            }
        }#}}}
        plot_antibody <- function(antibody) {#{{{
            print(RidgePlot(dt, features = antibody))
            print(VlnPlot(dt, features = antibody, pt.size = 0) +geom_boxplot(width=0.1,fill="white"))
        }#}}}
        for(feature in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"))
            ploti(feature)
        for(antibody in antibodies)
            plot_antibody(antibody)
    }#}}}
    qc_plot <- function(dt, dt_name, lim=NULL, suffix, traits, antibodies) {#{{{
        if(is.null(antibodies))
            antibodies <- grep('totalc', rownames(dt), value = TRUE, ignore.case = TRUE)
        pdf(paste0("results/0.QC_plots.", dt_name, suffix, ".pdf"), width=5, height=5)
        for(ident in c(dt@project.name, traits)) {
            Idents(dt) <- ident
            qc_plots0(dt, antibodies, lim)
        }
        dev.off()
    }#}}}
    qc_filter <- function(dt, lim) {#{{{
        my_filter <- function(val, lim) {#{{{
            if(length(lim) > 0) {
                res <- if(is.na(lim)) rep(TRUE, length(val)) else
                    val >= lim[1]
            }
            if(length(lim) > 1) res <- res & val <= lim[2]
            return(res)
        }#}}}
        if(length(lim) < 1) return(dt)
        keep <- mapply(my_filter, lapply(names(lim), function(x) dt[[x]]), lim)
        if(is.matrix(keep)) keep <- apply(keep, 1, all)
        return(dt[, keep])
    }#}}}
    qc_normalize <- function(dt, method = "sct") {#{{{
        imethod <- pmatch(method, c("sctransform", "scale"))
        if(imethod %in% 1)  {# SCT {{{
            dt <- SCTransform(dt, vars.to.regress = vars.to.regress, verbose = FALSE, ...)
        }#}}}
        if(imethod %in% 2) {# SCALE & NORMALIZE {{{
            dt <- dt %>%
                NormalizeData %>%
                FindVariableFeatures(selection.method = "vst", ...) %>%
                ScaleData(vars.to.regress = vars.to.regress)
        }#}}}
        saveRDS(dt, paste0("results/1.data_normalized.", dt_name, suffix, ".rds"))
        return(dt)
    }#}}}

    if(!test.na(percent_pattern)) dt <- calc_percentage(dt, percent_pattern)
    if(!test.na(qc_lim)) dt <- qc_filter(dt, qc_lim)
    if(!test.na(plot_lim)) qc_plot(dt, dt_name, plot_lim, suffix, traits, antibodies)
    if(!test.na(norm_method)) dt <- qc_normalize(dt, method = norm_method)

    return(dt)
}#}}}
my_DimReduct <- function(dt, dt_name, traits=c(), antibodies = NULL, dr = TRUE, suffix = "") {#{{{
    if(dr) {#{{{
        dt <- dt %>% RunPCA
        dt <- dt %>% RunUMAP(dims = 1:30)
        dt <- dt %>% FindNeighbors(dims = 1:30)
        dt <- dt %>% FindClusters(dims = 1:30)

        cat("Saving results after finding clusters ...")
        saveRDS(dt, paste0("results/2.clusters.", dt_name, suffix, ".rds"))
        cat("\n")
    } #}}}

    umap0 <- function(dt, antibodies, split.by = NULL) { #  {{{
        plot_antibody <- function(antibody) {#{{{
            print(VlnPlot(dt, features = antibody, split.by = split.by,
                          pt.size = 0, combine = FALSE) +
                  geom_boxplot(width=0.1,fill="white"))
        }#}}}
        print(DimPlot(dt, label = TRUE, reduction = "umap", pt.size = 0.3,
                      repel = TRUE, split.by = split.by, combine = FALSE))
        for(antibody in antibodies)
            plot_antibody(antibody)
    }#}}}
    umap_plot <- function(dt, dt_name, suffix, traits, antibodies) {#{{{
        if(is.null(antibodies))
            antibodies <- grep('totalc', rownames(dt), value = TRUE, ignore.case = TRUE)
        Idents(dt) <- "seurat_clusters"

        pdf(paste0("results/2.umap.", dt_name, suffix, ".pdf"), width=5, height=5)
        umap0(dt, antibodies, split.by = NULL)
        for(trait in traits)
            umap0(dt, antibodies, split.by = trait)
        dev.off()
    }#}}}

    umap_plot(dt, dt_name, suffix, traits, antibodies)
    return(invisible(dt))
}#}}}
my_DE <- function(dt, dt_name, traits) {#{{{
    DE<- function(dt, trait, name) {#{{{
        Idents(dt) <- trait
        cluster.markers <- FindAllMarkers(object = dt,
                                          only.pos = TRUE,
                                          min.pct = 0.10,
                                          return.thresh = 0.01,
                                          slot = 'data',
                                          logfc.threshold = 0.25)
        my.export.findmarkers(cluster.markers, paste0("results/2.DE_cluster_makers.", dt_name, ".", name, ".tsv"))
    }#}}}
    DE(dt, "seurat_clusters", "seurat_clusters")
    for(cluster in sort(unique(dt$seurat_clusters))) {#{{{
        dt_i <- dt[, dt$seurat_clusters %in% cluster]
        for(trait in traits) {#{{{
            DE(dt_i, trait, paste0(trait, ".within_cluster", cluster))
        }#}}}
    }#}}}

}#}}}
