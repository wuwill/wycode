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
setwd("/scratch/gtac/analysis/rna_seq/8998_3_magee-s4483")
rpkm <- get.gene_count('s4483_magee_unknown_mouse')
write.table(rpkm, file="all.gene_RPKM.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

