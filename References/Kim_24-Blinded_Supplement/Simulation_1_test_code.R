#' This file is used for the Section 6.1 Probability Proportional to Size Sampling with Replacement.

rm(list = ls())
library(rootSolve)
library(foreach)
library(parallel)
library(doParallel)

source("Simulation_1_functions.R")

# Test the code.
n = 50
N = 1000
beta1 = 1
proposed.bootstrap(seed.i = 1, N = N, n = n, beta1 = beta1)


# Do parallel
dir.create("Auxiliary/")
for (beta1 in c(1, 1.05, 1.1)) {
    for (n in c(75, 50)) {
        if (n == 50) {
            N = 1000
        } else if (n == 75) {
            N = 2000
        }
        print(paste(beta1, n, N))
        
        
        cl = makeCluster(10)
        registerDoParallel(cl)
        result.set = foreach(i = 1:1000, .combine = rbind, .packages = c("rootSolve"), .export = c("gen.popu", "sample.f", "estimation.f", "log.like.f", 
            "X.S.2.f"), .errorhandling = "remove") %dopar% {
            proposed.bootstrap(seed.i = i, N = N, n = n, beta1 = beta1)
        }
        stopCluster(cl)
        write.csv(result.set, file = paste0("Auxiliary/Simulation_1_N_", N, "_beta1_", 10 * beta1, ".csv"), row.names = FALSE)
    }
    
}


# Analyze the results The output corresponds to NLR, NQS, LS, BLR, and BQS, respectively.
for (beta1 in c(1, 1.05, 1.1)) {
    for (n in c(75, 50)) {
        if (n == 50) {
            N = 1000
        } else if (n == 75) {
            N = 2000
        }
        print(paste(beta1, n, N))
        result.set = read.csv(file = paste0("Auxiliary/Simulation_1_N_", N, "_beta1_", 10 * beta1, ".csv"))
        
        print(apply(result.set, 2, mean))
        
    }
    
}

