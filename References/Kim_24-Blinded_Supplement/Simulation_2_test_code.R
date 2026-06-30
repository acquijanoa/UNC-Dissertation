#' This simulation is conducted to test the performance of the bootstrap quasi-score method with the naive one in Section 6.2 Stratified Two-Stage Cluster Sampling

rm(list = ls())
library(plyr)
library(foreach)
library(parallel)
library(doParallel)
library(survey)

source("Simulation_2_functions.R")

p=c(1/4,1/8,1/8,1/8,1/16,1/16,1/8,1/16,1/16) # independent case
N=2000
n = 75
alpha = 0.05
one.iter.f(seed.i = 1, n = n, N = N, p = p, alpha = 0.05)

for (case in c(1,2,3))
{
    if (case == 1){ 
        p =c(1/4,1/8,1/8,1/8,1/16,1/16,1/8,1/16,1/16) # This is the  independent case
    } else if (case == 2){
        p=c(1/4,1.4/8,1.4/8,0.6/8,1/16,1.4/16,0.6/8,0.6/16,1/16) # Gamma = 0.0166
    } else if (case == 3){
        p=c(1/6,1/12,1/12,1/12,1/12,1/6,1/12,1/6,1/12) #Gamma = 0.125
    }
    for (N in c(2000,10000)){
        if (N == 10000){
            n=150
        } else if (N==2000){
            n=75
        }
        
        cl=makeCluster(10)
        registerDoParallel(cl)
        result=foreach(i=1:1000,.combine=rbind,.packages=c("survey"))%dopar%
            {
                one.iter.f(seed.i = i, n = n, N = N, p = p, alpha = 0.05)
            }
        stopCluster(cl)
        
        cat(N,n,case,":",round(apply(result,2,mean,na.rm=TRUE),2))
        print(" ")
    }
}
