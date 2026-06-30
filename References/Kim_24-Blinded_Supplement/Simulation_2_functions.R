#' This simulation is conducted to test independence in a two-way table

#' The function below is used to generate the population
gen.popu = function(N = 10000,p){
    #' N: population size
    #' p: probability vector
    
    popu = t(rmultinom(N,1,p))
    size = matrix(rexp(length(p),rate = 1)+0.5,ncol = 1)
    y = popu%*% size
    prob = abs(y)/sum(abs(y))
    popu = cbind(prob,popu)
    
    return(popu)
}

#' The function below is used to obtain the sample.
sample.f =  function(popu,n){
    #' popu: the finite population
    #' n: the sample size
    first.selected = sample(1:nrow(popu),size = n,prob = popu[,1],replace = TRUE) # Use PPS
    sample.result = popu[first.selected,]
    
    return(sample.result)
}


# The following is used to calculate the level of non-independence as well as the chi-squared test statistic
gamma.gen = function(p)
{
    #' p: the (estimated) probability vector
    
    if (round(sum(p),10) != 1)
        stop("Your probability vector is not normalized.")
    p.mat = matrix(p,3,3,byrow = TRUE)
    
    gamma.re = 0
    for (i in 1:3){
        for (j in 1:3){
            gamma.re = gamma.re + (p.mat[i,j] - sum(p.mat[i,]) * sum(p.mat[,j]))^2/(sum(p.mat[i,]) * sum(p.mat[,j]))
        }
    }
    return(gamma.re)
}

# The following is the function used to estimate the probability vector under PPS sampling
estimation.f = function(sample.result){
    #' sample.result: a sample with the corresponding selection probability in the first column
    
    weight.s = 1/sample.result[,1]
    data.s = sample.result[,-1]
    
    est.s = (matrix(weight.s,nrow = 1) %*% data.s)/nrow(sample.result)
    result.s = est.s/sum(est.s)
    
    return( result.s)
}

# The following is the function used to calculate the test statistic with respect to the second-order Rao-Scott correstion method
Rao_Scott_second_order = function(sample.result){
    #' sample.result: a sample with the corresponding selection probability in the first column
    
    row.colume = expand.grid(col = 1:3,row = 1:3)
    data.result = data.frame(w = 1/(sample.result[,1] * nrow(sample.result)),
                             col = NA,
                             row = NA)
    index = sample.result[,-1] %*% matrix(1:9,ncol = 1)
    data.result$col = row.colume$col[index]
    data.result$row = row.colume$row[index]
    data.result$fpc = N
    des = svydesign(id = ~1, weights = ~w,data = data.result, fpc =~fpc)
    aa = svychisq(~col + row,des)
    return(aa$p.value)
}



# The following is for one Monte Carlo simulation
one.iter.f = function(seed.i,n,N,p, alpha){
    #' seed.i: used to set seed.
    #' n: sample size
    #' N: population size
    #' p: probability vector 
    #' alpha: nominal significance level.
    
    set.seed(12356328+seed.i)
    result.compare = rep(NA,3) # This is used to store the simulation results.
    
    popu = gen.popu(N,p)
    sample.result = sample.f(popu = popu,n = n)
    
    estimates.ii = try(estimation.f(sample.result = sample.result),silent = TRUE)
    
    if (length(estimates.ii) >1){
        X.s.2.ori = gamma.gen(estimates.ii) * n
        
        if (is.double(X.s.2.ori)){
            result.compare[1] = X.s.2.ori > qchisq(1-alpha,df = 4)
            result.compare[2] = Rao_Scott_second_order(sample.result = sample.result)<alpha

            one.iter = function(i){ # Generate bootstrap sample and the bootstrap test statistic
                sample.result.jj = sample.result
                k.n = n/(n-1)
                wei.jj = k.n * rmultinom(1,n-1,rep(1/n,n))
                sample.result.jj[,1]= sample.result.jj[,1]/wei.jj
                
                
                if (! is.null(sample.result.jj)){
                    estimates.jj = estimation.f(sample.result = sample.result.jj)
                    
                    X.2.gen = function(estimates.ii,estimates.jj){
                        p.mat.ii = matrix(estimates.ii,3,3,byrow = TRUE)
                        p.mat.jj = matrix(estimates.jj,3,3,byrow = TRUE)
                        
                        gamma.re = 0
                        for (i in 1:3){
                            for (j in 1:3){
                                gamma.re = gamma.re + (p.mat.jj[i,j] - sum(p.mat.jj[i,]) * sum(p.mat.jj[,j]) - (p.mat.ii[i,j] - sum(p.mat.ii[i,]) * sum(p.mat.ii[,j])))^2/(sum(p.mat.ii[i,]) * sum(p.mat.ii[,j]))
                            }
                        }
                        return(gamma.re)
                    }
                    X.s.2.jj =  X.2.gen(estimates.ii,estimates.jj) * nrow(sample.result.jj)
                    
                    return(X.s.2.jj)
                } else{ 
                    stop("Error")}
            }
            
            
            result.i = rep(NA,2000)
            for (jjj in 1:2000){
                aa = try(one.iter(jjj),silent = TRUE)
                if (is.numeric(aa))
                    result.i[jjj]=aa
            }
            result.compare[3] = X.s.2.ori> quantile(result.i,1-alpha,na.rm = TRUE)
            
        }
    }
        
    
    return(result.compare)
    
}



