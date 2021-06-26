# Nuti G. An Efficient Algorithm for Bayesian Nearest Neighbours. Methodology and Computing in Applied Probability. 2019 Dec;21(4):1251-8.
bnn <- function(nn_class, alpha=10, beta=10, p_gamma = 0.05, class1=nn_class[1]) { # Bayesian Nearest Neighours {{{

    # nn_class: classes for nearest neighours in the order from nearest, to next nearest, ....
    #           By default, include 5 * (alpha + beta) neighours;
    # alpha, beta: Beta prior proablity parameters;
    # p_gamma: proablity of streak extension, ie, not a breakpoint at the next position
    # output: a list containing changepoint position, and the posterior proablity

    classes <- unique(nn_class)
    x <- rev(nn_class) # reorder: farest -> nearest
    x <- 2 - x %in% class1 # class1 -> 1; class2 -> 2
    x <- c(x, 1) # calculate proablity of observing class1 at after finishing updating
    n <- length(x)

    get.prop <- function(x) x / sum(x)
    eta <- matrix(NA, ncol = 2, nrow = n + 1) # beta coefficients for different streak lengths
    eta[1, ] <- eta0 <- c(alpha, beta) # priors

    prop0 <- get.prop(eta0)
    get.pi0 <- function(t) prop0[x[t]] # pi when streak length == 0

    pk <- matrix(0, nrow = n+1, ncol=n+1) # p(k_t = i | x)
    pkx <- matrix(0, nrow = n+1, ncol=n+1) # p(k_t = i,  x)
    pi <- matrix(0, nrow = n+1, ncol=n+1) # p(x_t | k_(t-1) = i,  x)
    pi[1, 1] <- pk[1, 1] <- pkx[1, 1] <- 1 # p(k0=0) <- 1

    for(t in 1:n){ # (t+1)-th row, (i+1)-th column of pk / pkx;  t starting from 1
        # observe the next x_t
        #   pi(i): p(x_t | k_(t-1) = i, eta_i)
        #   pk(i): p(k_t = k_(t-1) + 1, x | k_(t-1) = i)
        for(i in 0:(t-1)){
            # Compute predictive prob:
            #   pi_i = p(x_t | k_(t-1), eta_i)
            pi[t + 1, i + 1] <- get.prop(eta[i + 1, ])[x[t]]

            # Compute growth prob:
            #   p(k_t = k_(t-1) + 1, x_(0:t)) = p(k_(t-1), x_(0:t-1)) * pi_t * p_gamma (ERROR?)
            pkx[t + 1, i + 2] <- pkx[t, i + 1] * pi[t+1, i + 1] * (1 - p_gamma)
        }
        if(t == n) {
            posterior = sum(pk[n, ] * pi[n + 1, ])
            return(list(pi = pi, pk = pk, posterior = posterior, bayes_factor = posterior / (1 - posterior), class1 = class1))
        }

        # Compute change-point prob:
        #   p(kt =0, x_0:t) = sum_{} 
        pkx[t + 1, 1] <- sum( pkx[t, 1:t]  * p_gamma) * pi[t + 1, 1]

        # Compute evidence
        #   p(x_0:t) = sum_kt { p(kt, x_0:t) }
        px <- sum(pkx[t + 1, 1:(t + 1)])

        eta.prev <- eta
        for(i in 0:t){
            # Compute prob of k:
            #   p(k_t = i | x_0:t) = p(k_t = i, x_0:t) / p(x_0:t)
            pk[t + 1, i + 1] = pkx[t + 1, i + 1] / px

            # Update distributions:
            #   eta_i <<- x_t
           if(i==0) eta[1, ] <- eta0 else {
               eta[i + 1, ] <- eta.prev[i, ]
               eta[i + 1, x[t]] <- eta[i + 1, x[t]] + 1
           }
        }
    }

    # return: p(k | x)
} #}}}
get_sample_weight <- function(x, group = NULL) {#{{{
    unique_x <- sort(setdiff(x, NA))
    res <- rep(1, length(x))
    res[is.na(x)] <- NA
    if(is.null(group)) group <- res
    if(length(group) == 1 && is.na(group)) return(res)

    groups <- setdiff(group, NA)
    for(g in groups) {#{{{
        i <- group %in% g
        g_count <- sum(i)
        nc <- length(setdiff(x[i], NA))
        for(v in unique_x) {#{{{
            j <- i & x %in% v
            if(any(j)) res[j] <- g_count / nc / sum(j)
        }#}}}
    }#}}}

    return(res)
}#}}}
bnn0 <- function(nn_class, alpha=10, beta=10, p_gamma = 0.05, class1=nn_class[1], weights = NULL, group = NULL) { # weight samples {{{
    # weights = 1 - unweighted
    #         = NULL - weight by sample frequency for unblanced data
    # nn_class: classes for nearest neighours in the order from nearest, to next nearest, ....
    #           By default, include 5 * (alpha + beta) neighours;
    # alpha, beta: Beta prior proablity parameters;
    # p_gamma: proablity of streak extension, ie, not a breakpoint at the next position
    # output: a list containing changepoint position, and the posterior proablity
    #       - pi, prob(x_t | k), proablity of the current position being class1 for different run lengths
    #       - pk, prob(k | x), proablity of being a run lengths
    #       - posterior, posterior proablity of the current posistion is from class1

    x0 <- rev(nn_class) # reorder: farest -> nearest
    classes <- c(class1, setdiff(x0, class1))
    x0 <- match(x0, classes)
    nc <- length(classes)
    if(is.null(weights)) {#{{{
        weights <-  get_sample_weight(x0, group = group)
    } else if(length(weights) == 1 && is.na(weights)) {
        weights <-  get_sample_weight(x0, group = NA)
    } else if(length(weights) < length(x0)) {
        weights <- rep(weights, length = x0)
    }#}}}

    x <- 2 - x0 %in% 1 # class1 -> 1; class2 -> 2
    x <- c(x, 1) # add 1 to the last so as to calculate proablity of observing class1 at after finishing updating
    n <- length(x)
    x0 <- c(x0, 1)

    eta <- eta0 <- rbind(c(alpha, beta)) # priors
    eta.sum <- eta.sum0 <- sum(eta0)
    pk.prev <- pkx.prev <-  1 # pk = prob(k | x); pkx = prob(k, x)

    for(t in 1:n){
        # observe the next x_t
        #   pi(i): p(x_t | k_(t-1) = i, eta_i)
        #   pk(i): p(k_t = k_(t-1) + 1, x | k_(t-1) = i)

        # Compute predictive prob:
        #   pi_i = p(x_t | k_(t-1), eta_i)
        pi <- eta[,x[t]] / eta.sum

        if(t == n) {
            posterior = sum(pk.prev * pi, na.rm=TRUE)
            return(list(pi = pi, pk = pk.prev, posterior = posterior, bayes_factor = posterior / (1 - posterior), class1 = class1))
        }

        # Compute change-point prob:
        #   p(kt =0, x_0:t) = sum_{} 
        # Compute growth prob:
        #   p(k_t = k_(t-1) + 1, x_(0:t)) = p(k_(t-1), x_(0:t-1)) * pi_t * p_gamma (ERROR?)
        pkx <- c(sum(pkx.prev)  * p_gamma * pi[1], pkx.prev * pi * (1 - p_gamma))

        # Compute evidence
        #   p(x_0:t) = sum_kt { p(kt, x_0:t) }
        px <- sum(pkx)
        # Compute prob of k:
        #   p(k_t = i | x_0:t) = p(k_t = i, x_0:t) / p(x_0:t)
        pk <- pkx / px

        # Update distributions:
        #   eta_i <<- x_t
        incr <- weights[x0[t]]
        eta[, x[t]] <- eta[, x[t]] + incr
        eta <- rbind(eta0, eta)
        eta.sum <- c(eta.sum0, eta.sum + incr)
        pk.prev <- pk
        pkx.prev <- pkx
    }

} #}}}
bnn2 <- function(nn_class, alpha=10, beta=10, p_gamma = 0.05, class1=nn_class[1]) { # resamples {{{
    # nn_class: classes for nearest neighours in the order from nearest, to next nearest, ....
    #           By default, include 5 * (alpha + beta) neighours;
    # alpha, beta: Beta prior proablity parameters;
    # p_gamma: proablity of streak extension, ie, not a breakpoint at the next position
    # output: a list containing changepoint position, and the posterior proablity
    #       - pi, prob(x_t | k), proablity of the current position being class1 for different run lengths
    #       - pk, prob(k | x), proablity of being a run lengths
    #       - posterior, posterior proablity of the current posistion is from class1

    classes <- unique(nn_class)
    x <- rev(nn_class) # reorder: farest -> nearest
    x <- 2 - x %in% class1 # class1 -> 1; class2 -> 2
    x <- c(x, 1) # add 1 to the last so as to calculate proablity of observing class1 at after finishing updating
    n <- length(x)

    eta <- eta0 <- rbind(c(alpha, beta)) # priors
    eta.sum <- eta.sum0 <- sum(eta0)
    pk.prev <- pkx.prev <-  1 # pk = prob(k | x); pkx = prob(k, x)

    for(t in 1:n){
        # observe the next x_t
        #   pi(i): p(x_t | k_(t-1) = i, eta_i)
        #   pk(i): p(k_t = k_(t-1) + 1, x | k_(t-1) = i)

        # Compute predictive prob:
        #   pi_i = p(x_t | k_(t-1), eta_i)
        pi <- eta[,x[t]] / eta.sum

        if(t == n) {
            posterior = sum(pk.prev * pi, na.rm=TRUE)
            return(list(pi = pi, pk = pk.prev, posterior = posterior, bayes_factor = posterior / (1 - posterior), class1 = class1))
        }

        # Compute change-point prob:
        #   p(kt =0, x_0:t) = sum_{} 
        # Compute growth prob:
        #   p(k_t = k_(t-1) + 1, x_(0:t)) = p(k_(t-1), x_(0:t-1)) * pi_t * p_gamma (ERROR?)
        pkx <- c(sum(pkx.prev)  * p_gamma * pi[1], pkx.prev * pi * (1 - p_gamma))

        # Compute evidence
        #   p(x_0:t) = sum_kt { p(kt, x_0:t) }
        px <- sum(pkx)
        # Compute prob of k:
        #   p(k_t = i | x_0:t) = p(k_t = i, x_0:t) / p(x_0:t)
        pk <- pkx / px

        # Update distributions:
        #   eta_i <<- x_t
        eta[, x[t]] <- eta[, x[t]] + 1
        eta <- rbind(eta0, eta)
        eta.sum <- c(eta.sum0, eta.sum + 1)
        pk.prev <- pk
        pkx.prev <- pkx
    }

} #}}}
get.local_bnn <- function(dat=NULL, bnn_input, chngpnt = NULL, n_ave = 5, top = TRUE, n_top = 1, name = "bnn_top", nn = NULL) {#{{{
    if(!top) bnn_input <- 1 - bnn_input
    if(is.null(nn)) {
        library(FNN)
        nn <- get.knn(dat, k = n_ave)#nn.index; nn.dist
    }
    n <- nrow(nn$nn.index)

    ave_fun <- function(x) exp(mean(log(x)))
    ave_bnn <- apply(matrix(bnn_input[cbind(1:n, nn$nn.index[, 1:n_ave])], nrow = n), 1, ave_fun)

    i <- NULL; d <- neighors <- as.data.frame(matrix(0, nrow = n, ncol = n_top))
    while(length(i) < n_top) {
        i_ <- which.max(ave_bnn)
        d_ <- sqrt(apply(sweep(dat, 2, dat[i_,], "-")^2, 1, sum))
        # masks neighors of i_
        ii <- order(d_) # index of neighors of i_
        j <- if(is.null(chngpnt)) {
            max(n_ave * 3, min(which(ave_bnn[ii] < median(ave_bnn, na.rm = TRUE)), na.rm = TRUE))
        } else chngpnt[i_] + 2 * n_ave
        ave_bnn[ii[1:j]] <- NA

        i <- c(i, i_)
        d[, length(i)] <- d_
        neighors[ii[1:(j - n_ave)], length(i)] <- 1
    }

    colnames(d) <- paste0(name, 1:n_top, "_d", i)
    colnames(neighors) <- paste0(name, 1:n_top, "_n", i)
    
    res <- if(is.null(name))
        list(i = i, d = d, n = neighors) else
            cbind(d, neighors)
    # i - which is selected as the top / bottom;
    # d - distance to the top / bottom
    # n - neighorhood of the top / bottom
    return(res)
}#}}}
test_bnn__ <- function() {#{{{
    set.seed(0)
    res <- list(x = list(),
                prob1 = list(),
                prob2 = list(),
                n1 = list(),
                n2 = list(),
                bnn = list()
    )
    simulate <- function(i, ...) {#{{{
        res <- res
        input <- list(...)
        assign <- function(term) {#{{{
            if(is.null(input[[term]]) && i > 1) {
                res[[term]][[i]] <<- res[[term]][[i - 1]]
            } else {
                res[[term]][[i]] <<- input[[term]]
            }
        }#}}}
        for(term in c("n1", "n2", "prob1", "prob2", "x")) assign(term)
        if(any(names(input) %in% c("n1", "n2"))) res$x[[i]] <- c(sample(2, res$n1[[i]], TRUE, res$prob1[[i]]),
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
        if(any(names(input) %in% c("prob1"))) res$x[[i]][1:res$n1[[i]]] <- sample(2, res$n1[[i]], TRUE, res$prob1[[i]])
        if(any(names(input) %in% c("prob2"))) res$x[[i]][res$n1[[i]] + 1:res$n2[[i]]] <- sample(2, res$n2[[i]], TRUE, res$prob2[[i]])
        return(res)
    }#}}}

    # null case with equal chance of two groups
    i <- 1
    res <- simulate(i, prob1 = c(0.5, 0.5), prob2 = c(0.1, 0.9), n1 = 10, n2 = 20)
    res$bnn[[i]] <- bnn0(res$x[[i]], 1, 1)

    # more chance of being 2
    i <- 2
    res <- simulate(i, prob1 = c(0.1, 0.9), prob2 = c(0.5, 0.5))
    res$bnn[[i]] <- bnn0(res$x[[i]], 1, 1)

    # more chance of being 2, with biased background
    i <- 3
    res <- simulate(i, prob2 = c(0.8, 0.2))
    res$bnn[[i]] <- bnn0(res$x[[i]], 1, 1)


    # more chance of being 2
    i <- 3
    res$prob1[[i]] <- res$prob1[[i-1]]
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn0(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn0(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn0(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn0(res$x[[i]])

    return(res)
}#}}}
