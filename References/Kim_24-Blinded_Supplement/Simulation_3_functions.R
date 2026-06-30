#' This simulation is conducted to test the performance of the bootstrap quasi-score method with the naive one in Section 6.2 Stratified Two-Stage Cluster Sampling

#' The function below is used to generate the population
gen.popu = function(N.strata = 50, N.1 = 30, N.2 = 10) {
    #' N.strata: number of stratum
    #' N.1: number of clusters in each stratum
    #' N.2: number of elements in each cluster
    
    stratum.effect = rexp(N.strata, rate = 1)
    N.cluster = rpois(N.strata, lambda = stratum.effect) * 5 + N.1  # This is the number of clusters in each stratum.
    
    popu = NULL
    for (s in 1:N.strata) {
        cluster.effect = rexp(N.cluster[s])
        N.ele = rpois(sum(N.cluster[s]), lambda = stratum.effect[s] + cluster.effect) * 5 + N.2
        
        popu.s = data.frame(stratum = s, cluster = rep(1:N.cluster[s], N.ele), x = runif(sum(N.ele), min = 0, max = 20), M.h = N.cluster[s], M.hi = rep(N.ele, 
            N.ele), N.h = sum(N.ele))  # stratum index, cluster index, covariate, number of clusters, cluster size, stratum size
        popu.s$y = popu.s$x + 0.3 * rep(cluster.effect, N.ele) + 0.3 * stratum.effect[s] + rnorm(nrow(popu.s))
        
        popu = rbind(popu, popu.s)
    }
    
    popu$N = nrow(popu)  # Population size
    
    return(popu)
}


#' The function below is used to obtain the sample.
sample.f = function(popu, n.1 = 5, n.2 = 5) {
    #' popu: finite population
    #' n.1: number of sampled clusters
    #' n.2: number of sampled units
    
    # Initialize the sample
    sample.result = popu[rep(1, max(popu$stratum) * n.1 * n.2), ] + NA
    sample.result$w1 = NA
    sample.result$w2 = NA
    sample.result$weight = NA
    sample.result$index = NA
    
    
    index.now = 1
    inc.step = n.2
    for (s in unique(popu$stratum)) {
        popu.s = subset(popu, stratum == s)
        summary.s = ddply(popu.s, "cluster", summarize, l = unique(M.hi))
        
        sampled.clus = sample(1:nrow(summary.s), size = n.1, prob = summary.s$l, replace = TRUE)
        for (i in 1:n.1) {
            sampled.i = subset(popu.s, cluster == sampled.clus[i])
            sampled.ele.i = sample(1:nrow(sampled.i), n.2, replace = FALSE)
            
            sample.result.i = sampled.i[sampled.ele.i, ]
            sample.result.i$w1 = (n.1 * unique(sample.result.i$M.hi)/unique(sample.result.i$N.h))^(-1)
            sample.result.i$w2 = unique(sample.result.i$M.hi)/n.2
            sample.result.i$weight = sample.result.i$w1 * sample.result.i$w2
            sample.result.i$index = i
            sample.result[index.now:(index.now + inc.step - 1), ] = sample.result.i
            index.now = index.now + inc.step
        }
    }
    
    return(sample.result)
}

#' The following is the estimating function of beta. 
estimation.f = function(sample.result, beta1) {
    #' sample.result: the generated sample
    #' beta1: the null hypothesis
    
    N = unique(sample.result$N)
    mat.x = model.matrix(y ~ x, data = sample.result)
    W = diag(sample.result$weight)
    
    # Obtain the pseudo maximum loglikelihood parameter.
    theta.full = solve(t(mat.x) %*% W %*% mat.x) %*% t(mat.x) %*% W %*% matrix(sample.result$y, ncol = 1)
    
    # Obtain the pseudo maximum loglikelihood parameter under the null hypothesis.
    parameter.solve.constraint = function(x) {
        beta.h = matrix(c(x[1], beta1), ncol = 1)
        resi.h = (sample.result$y - mat.x %*% beta.h)
        c(F1 = sum(sample.result$weight * resi.h)/N)
    }
    b1 = sum((sample.result$y - beta1 * sample.result$x) * sample.result$weight)/sum(sample.result$weight)
    theta.constraint = c(b1, beta1)
    
    return(list(theta.full = theta.full, theta.cons = theta.constraint))
}


#' The following is used to calculate the quasi score.
X.S.2.f = function(sample.result, theta) {
    #' sample.result: the generated sample
    #' beta: the estimated parameter under the null hypothesis
    
    N = unique(sample.result$N)
    
    X.design.matrix = t(model.matrix(~sample.result$x))
    resi.x = sample.result$y - t(X.design.matrix) %*% matrix(theta[1:2], ncol = 1)
    W = diag(sample.result$weight)
    
    # Calculate the pseudo Fisher information.
    I = matrix(0, 2, 2)
    I[1:2, 1:2] = X.design.matrix %*% W %*% t(X.design.matrix)
    I = I/N
    I.w2.1 = I[2, 2] - matrix(I[2, -2], nrow = 1) %*% solve(I[-2, -2]) %*% matrix(I[2, -2], ncol = 1)
    S.w2.0 = sum(sample.result$weight * resi.x * sample.result$x)/N
    
    return(S.w2.0 * I.w2.1^(-1) * S.w2.0)
}



proposed.bootstrap = function(seed.i, N.strata, N.1, N.2, n.1, n.2, beta1, alpha = 0.05) {
    #' seed.i: used to set seed.
    #' N.strata: number of stratum
    #' N.1: number of clusters in each stratum
    #' N.2: number of elements in each cluster
    #' n.1: number of sampled clusters
    #' n.2: number of sampled units
    #' beta1: null hypothesis
    #' alpha: nominal significance level.
    
    # Generation of the finite population and sample.
    set.seed(1234 + seed.i)
    popu = gen.popu(N.strata = N.strata, N.1 = N.1, N.2 = N.2)
    sample.result.ii = sample.f(popu = popu, n.1 = n.1, n.2 = n.2)
    
    # The pseudo maximum likelihood estimator and the one under null hypothesis
    estimates.ii = estimation.f(sample.result = sample.result.ii, beta1 = beta1)
    theta.full.ii = estimates.ii$theta.full
    theta.cons.ii = estimates.ii$theta.cons
    
    # Obtain the test statistics
    X.s.2.ori = X.S.2.f(sample.result = sample.result.ii, theta = theta.cons.ii)
    
    result.all = rep(NA, 2)  # This variable is used to store the testing result using the two bootstrap methods.
    result.all[1] = nrow(sample.result.ii) * X.s.2.ori > qchisq(1 - alpha, df = 1)  # Naive Quasi-score method
    
    
    one.iter = function(seed.j) {
        #' seed.j: used to set seed
        
        set.seed(1341232 + seed.i * 10 + seed.j)
        
        S.j = sample.result.ii  # Generate the bootstrap weights.
        k.h = n.1/(n.1 - 1)  # The weight adjustment factor.
        reweight.j = rep(NA, nrow(S.j))
        
        index.j = 1
        inc.step = n.1 * n.2
        for (iii in unique(S.j$stratum)) {
            reweight.j[index.j:(index.j + inc.step - 1)] = rep(rmultinom(1, n.1 - 1, rep(1/n.1, n.1)), each = n.2)
            index.j = index.j + inc.step
        }
        S.j$weight = S.j$weight * k.h * reweight.j
        
        estimates.jj = estimation.f(sample.result = S.j, beta1 = theta.full.ii[2])
        theta.full.jj = estimates.jj$theta.full
        theta.cons.jj = estimates.jj$theta.cons
        
        X.s.2.jj = X.S.2.f(sample.result = S.j, theta = theta.cons.jj)
        
        return(c(X.s.2.jj[1]))
    }
    
    result.i = rep(NA, 1000)
    for (seed.j in 1:1000) {
        result.i[seed.j] = one.iter(seed.j)
    }
    result.all[2] = X.s.2.ori[1] > quantile(result.i, 1 - alpha, na.rm = TRUE)
    
    return(result.all)  # The output corresponds to NQS and BQS, respectively
}
