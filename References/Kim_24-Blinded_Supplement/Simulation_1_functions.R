#' This file is used for the Section 6.1 Probability Proportional to Size Sampling with Replacement.

#' The function below is used to generate the population
gen.popu = function(N) {
    #' N: population size

    popu = data.frame(x = runif(N, -5, 5))
    mean_y = 1 +  popu$x
    popu$y =  mean_y + rnorm(N,mean=0,sd = 2) 
    ind_error_y_less1 = (popu$y - mean_y)*(popu$x)>0
    cov_xy = 6*ind_error_y_less1+1  
    popu$sel.prob = cov_xy 
    popu$sel.prob = popu$sel.prob/sum(popu$sel.prob)
    
    return(popu)
}

#' The function below is used to obtain the sample.
sample.f = function(popu, n) {
    #' popu: finite population
    #' n: sample size
    
    first.selected = sample(1:nrow(popu), size = n, prob = popu$sel.prob, replace = TRUE)  # Use PPS
    sample.result = popu[first.selected, ]
    sample.result$N = nrow(popu)
    sample.result$weight = 1/(n * sample.result$sel.prob)
    
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
    beta.hat = solve(t(mat.x) %*% W %*% mat.x) %*% t(mat.x) %*% W %*% matrix(sample.result$y, ncol = 1)
    sigma2.hat = sum(sample.result$weight)^(-1) * t(matrix(sample.result$y, ncol = 1) - mat.x %*% beta.hat) %*% W %*% (matrix(sample.result$y, ncol = 1) - 
        mat.x %*% beta.hat)
    theta.full = c(beta.hat, sigma2.hat)
    
    # Obtain the pseudo maximum loglikelihood parameter under the null hypothesis.
    parameter.solve.constraint = function(x) {
        beta.h = matrix(c(x[1], beta1), ncol = 1)
        resi.h = (sample.result$y - mat.x %*% beta.h)
        c(F1 = sum(sample.result$weight * resi.h)/N, F2 = sum(sample.result$weight * (-x[2] + resi.h^2))/N)
    }
    b1 = multiroot(parameter.solve.constraint, start = c(beta.hat[1], sigma2.hat))
    theta.constraint = c(b1$root[1], beta1, b1$root[2])
    
    return(list(theta.full = theta.full, theta.cons = theta.constraint))
}

#' The following function is used to calculate the weighted log-likelihood.
log.like.f = function(sample.result, theta) {
    #' sample.result: the generated sample
    #' beta: the estimated parameter
    
    N = unique(sample.result$N)
    
    pll = -1/2 * log(theta[3]) - (sample.result$y - theta[1] - sample.result$x * theta[2])^2/2/theta[3]
    
    return(sum(sample.result$weight * pll)/N)
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
    I = matrix(0, 3, 3)
    I[1:2, 1:2] = X.design.matrix %*% W %*% t(X.design.matrix)/theta[3]
    I[1:2, 3] <- I[3, 1:2] <- 1/theta[3]^2 * X.design.matrix %*% W %*% resi.x
    I[3, 3] = -sum(sample.result$weight * (1/2/theta[3]^2 - resi.x^2/theta[3]^3))
    I = I/N
    I.w2.1 = I[2, 2] - matrix(I[2, -2], nrow = 1) %*% solve(I[-2, -2]) %*% matrix(I[2, -2], ncol = 1)
    S.w2.0 = sum(sample.result$weight * resi.x * sample.result$x/theta[3])/N
    
    return(c(nrow(sample.result) * S.w2.0 * I.w2.1^(-1) * S.w2.0))  # Return the quasi score test statistic 
}

#' The following is used to implement the Lumley-Scott method.
L.S.method = function(sample.result, theta) {
    #' sample.result: the generated sample
    #' beta: the estimated parameter
    
    N = unique(sample.result$N)
    
    X.design.matrix = t(model.matrix(~sample.result$x))
    resi.x = sample.result$y - t(X.design.matrix) %*% matrix(theta[1:2], ncol = 1)
    W = diag(sample.result$weight)
    I = matrix(0, 3, 3)
    I[1:2, 1:2] = X.design.matrix %*% W %*% t(X.design.matrix)/theta[3]
    I[1:2, 3] <- I[3, 1:2] <- 1/theta[3]^2 * X.design.matrix %*% W %*% resi.x
    
    I[3, 3] = -sum(sample.result$weight * (1/2/theta[3]^2 - resi.x^2/theta[3]^3))
    I = I/N
    I.w2.1 = I[2, 2] - matrix(I[2, -2], nrow = 1) %*% solve(I[-2, -2]) %*% matrix(I[2, -2], ncol = 1)
    
    S.hat.theta = matrix(NA, nrow = nrow(sample.result), ncol = 3)
    for (i in 1:nrow(sample.result)) {
        S.hat.theta[i, ] = 1/N * sample.result$weight[i] * c(resi.x[i]/theta[3], resi.x[i] * sample.result$x[i]/theta[3], resi.x[i]^2/2/theta[3]^2 - 1/2/theta[3])
    }
    
    
    design.effect = (var(S.hat.theta) * nrow(sample.result) * nrow(sample.result))[2, 2]/(var(N * S.hat.theta/sample.result$weight) * (1 - nrow(sample.result)/N))[2, 
        2]
    var.theta = (solve(I) %*% var(S.hat.theta) %*% t(solve(I))) * nrow(sample.result) * nrow(sample.result)
    
    delta.hat = I.w2.1 * var.theta[2, 2]
    
    
    return(c(delta.hat, nrow(sample.result)/design.effect))  # Return the estimated delta and effective sample size.
}

#' The following function is used to do the simulation 1.
proposed.bootstrap = function(seed.i, N, n, beta1, alpha = 0.05) {
    #' seed.i: used to set seed.
    #' N: population size
    #' n: sample size
    #' beta1: null hypothesis
    #' alpha: nominal significance level.
    
    # Generation of the finite population and sample.
    set.seed(seed.i + 1235)
    popu = gen.popu(N)
    sample.result.ori = sample.f(popu = popu, n = n)
    
    # The pseudo maximum likelihood estimator and the one under null hypothesis
    estimates.ori = estimation.f(sample.result = sample.result.ori, beta1 = beta1)
    
    theta.full.ori = estimates.ori$theta.full  # The pseudo maximum likelihood estimator
    theta.cons.ori = estimates.ori$theta.cons  # The pseudo maximum likelihood estimator under the null hypothesis
    
    like.full.ori = log.like.f(sample.result = sample.result.ori, theta = theta.full.ori)  # The pseudo loglikelihood using the pseudo maximum likelihood estimator
    like.cons.ori = log.like.f(sample.result = sample.result.ori, theta = theta.cons.ori)  # The pseudo loglikelihood using the pseudo maximum likelihood estimator under the null hypothesis
    
    # Obtain the test statistics
    W.n.ori = -2 * nrow(sample.result.ori) * (like.cons.ori - like.full.ori)  # Wald statistic
    X.s.2.ori = X.S.2.f(sample.result = sample.result.ori, theta = theta.cons.ori)  # Quasi-score statistic
    delta.hat = L.S.method(sample.result = sample.result.ori, theta = theta.full.ori)  # Quantities used for Lumley-Scott method
    test.stat = W.n.ori/delta.hat[1]
    
    result.all = rep(NA, 5)  # This variable is used to store the testing result using the two bootstrap methods.
    
    result.all[1] = W.n.ori > qchisq(1 - alpha, df = 1)  # Naive Wald method
    result.all[2] = X.s.2.ori[1] > qchisq(1 - alpha, df = 1)  # Naive Quasi-score method
    result.all[3] = test.stat > qf(1 - alpha, df1 = 1, df2 = delta.hat[2] - 3)  # Lumley-Scott method
    
    
    one.iter = function(seed.j) {
        #' seed.j: used to set seed
        
        set.seed(1341232 + seed.i * 10 + seed.j)
        
        sample.result.boo = sample.result.ori  # Copy the original sample in order to change the bootstrap weights
        k.n = n/(n - 1)  # The weight adjustment factor.
        wei.boo = k.n * rmultinom(1, n - 1, rep(1/n, n))  # Use a multinomial distribution to generate bootstrap weights.
        sample.result.boo$weight = sample.result.boo$weight * c(wei.boo)  # Generate bootstrap weights.
        
        estimates.boo = estimation.f(sample.result = sample.result.boo, beta1 = theta.full.ori[2])
        theta.full.boo = estimates.boo$theta.full
        theta.cons.boo = estimates.boo$theta.cons
        
        
        like.full.boo = log.like.f(sample.result = sample.result.boo, theta = theta.full.boo)
        like.cons.boo = log.like.f(sample.result = sample.result.boo, theta = theta.cons.boo)
        
        W.n.boo = -2 * nrow(sample.result.boo) * (like.cons.boo - like.full.boo)
        X.s.2.boo = X.S.2.f(sample.result = sample.result.boo, theta = theta.cons.boo)
        
        return(c(W.n.boo, X.s.2.boo))
    }
    
    
    result.i = matrix(NA, 1000, 2)
    for (seed.j in 1:1000) {
        result.i[seed.j, ] = one.iter(seed.j)
    }
    
    result.all[4] = W.n.ori > quantile(result.i[, 1], 1 - alpha, na.rm = TRUE)  # Bootstrap Wald method
    result.all[5] = X.s.2.ori > quantile(result.i[, 2], 1 - alpha, na.rm = TRUE)  # Bootstrap quasi-score method
    
    return(result.all)  # The output corresponds to NLR, NQS, LS, BLR, and BQS, respectively.
}
