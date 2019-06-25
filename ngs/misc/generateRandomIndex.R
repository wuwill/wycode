generateRandomIndex <- function(n, nchar=ceiling(log(n, 4)*2), prev=NULL){ #{{{
    randomIndex <- if(is.null(prev)) c() else prev
    while(length(randomIndex)<n){ #{{{
        randomSeq <- matrix(sample(c('A', 'C', 'G', 'T'), n * nchar * nchar, replace=TRUE), ncol=nchar)
        randomIndex <- c(prev, setdiff(c(randomIndex, apply(randomSeq, 1, paste, collapse='')), prev))
    } #}}}
    return(randomIndex[1:n])
}#}}}
