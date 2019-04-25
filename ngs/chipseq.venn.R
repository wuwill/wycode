library(ChIPpeakAnno)
library(eulerr)
library(utils)
plot.eulerr <- function(venn){ #{{{
    nc <- ncol(venn)
    nn <- colnames(venn)[-nc]
    nv <- nc - 1

    ret <- list()
    #for(n in 1:nv){ #{{{
        #all.comb <- combn(1:nv, n)
        #for(i in 1:ncol(all.comb)){ #{{{
            #ii <- all.comb[,i]
            #ret[[paste0(nn[ii], collapse="&")]] <- sum(venn[apply(venn[,ii,drop=FALSE]==1, 1, all ),nc])
        #} #}}}
    #} #}}}
    for(i in 2:nrow(venn)){ #{{{
        ret[[paste(nn[venn[i, -nc] %in% 1], collapse="&")]] <- venn[i, nc]
    } #}}}
    counts <- unlist(ret)
    names(counts) <- names(ret)
    ret <- euler(counts)
    plot(ret, quantities=TRUE)
    return(ret)
} #}}}
my.venn.example <- function(peak_list){ #{{{
    v1 <- makeVennDiagram(peak_list, NameOfPeaks=names(peak_list), main="")
    p1 <- plot.eulerr(v1$vennCounts)
    print(plot(p1, quantities=TRUE, main=""))
} #}}}
