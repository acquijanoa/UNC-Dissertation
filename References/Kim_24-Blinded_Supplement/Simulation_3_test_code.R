#' This simulation is conducted to test the performance of the bootstrap quasi-score method with the naive one in Section 6.2 Stratified Two-Stage Cluster Sampling

rm(list = ls())
library(plyr)
library(rootSolve)
library(foreach)
library(parallel)
library(doParallel)

source("Simulation_2_functions.R")


# Test the code.
N.strata = 20
N.1 = 20
N.2 = 30
n.1 = 5
n.2 = 10
beta1 = 1
proposed.bootstrap(seed.i = 1, N.strata = N.strata, N.1 = N.1, N.2 = N.2, n.1 = n.1, n.2 = n.2, beta1 = beta1)

# Do parallel
for (beta1 in c(1, 1.01)) {
    for (N.strata in c(50, 20)) {
        N.1 = 20
        N.2 = 30
        n.1 = 5
        n.2 = 10
        print(paste(beta1, N.strata))
        
        cl = makeCluster(20)
        registerDoParallel(cl)
        result.set = foreach(i = 1:1000, .combine = rbind, .packages = c("plyr", "rootSolve"), .export = c("gen.popu", "sample.f", "estimation.f", "X.S.2.f"), 
            .errorhandling = "remove") %dopar% {
            proposed.bootstrap(seed.i = i, N.strata = N.strata, N.1 = N.1, N.2 = N.2, n.1 = n.1, n.2 = n.2, beta1 = beta1)
        }
        
        write.csv(result.set, file = paste0("Auxiliary/Simulation_2_N_strata_", N.strata, "_beta1_", 100 * beta1, ".csv"), row.names = FALSE)
    }
}


# The output corresponds to NQS and BQS, respectively
for (beta1 in c(1, 1.01)) {
    for (N.strata in c(50, 20)) {
        N.1 = 20
        N.2 = 30
        n.1 = 5
        n.2 = 10
        print(paste(beta1, N.strata))
        
        
        result.set = read.csv(file = paste0("Auxiliary/Simulation_2_N_strata_", N.strata, "_beta1_", 100 * beta1, ".csv"))
        print(apply(result.set, 2, mean))
    }
}
