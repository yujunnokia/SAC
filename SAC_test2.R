# Test the finite mixture of SACs model with multiple features
# 
# Author: Jun Yu
# Version: Oct 2013
###############################################################################


rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/SAC")
setwd("/nfs/guille/tgd/wonglab/yuju/SAC")
source("SAC.R")

# generate synthetic data
K <- 3 # number of groups
nBirder <- 1000 # number of birders
nCL <- 100 # number of checklists per birder
nFeatures <- 2

sigma <- 1
pi <- c(0.3, 0.3, 0.4)
beta <- matrix(c(2.5,1.5,0.5,1.0,1.0,0.5,-0.5,0.5,0.5),nrow=K,byrow=TRUE)

X <- matrix(sqrt(runif(nBirder*nCL*nFeatures, min=0.001, max=100)),ncol=nFeatures)
Y <- array(0,nBirder*nCL)
index <- array(0,nBirder*nCL)
for (i in 1:nBirder) {
    k <- sample(1:K,1,prob=pi)
    
    for (j in 1:nCL) {
        Y[(i-1)*nCL+j] <- rnorm(1,c(1,X[(i-1)*nCL+j,]) %*% beta[k,],sigma)
        index[(i-1)*nCL+j] <- i
    }
}
B <- list()
for (i in 1:nBirder) {
    B[[i]] <- which(index == i)
}

# learning
lambda <- 1
model <- RandomRestartEM.SAC(K,X,Y,B,lambda,2) 
cat("Output true model:\n")
print(beta)
print(pi)
print(sigma)




