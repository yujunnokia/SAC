# Test the finite mixture of SACs model with single feature
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

sigma <- 1
pi <- c(0.3, 0.3, 0.4)
beta <- matrix(c(2.5,1.5,1.0,1.0,-0.5,0.5),nrow=K,byrow=TRUE)

# generate training data
X <- matrix(sqrt(runif(nBirder*nCL, min=0.001, max=100)),ncol=1)
Y <- array(0,nBirder*nCL)
index <- array(0,nBirder*nCL)
for (i in 1:nBirder) {
    k <- sample(1:K,1,prob=pi)
    
    for (j in 1:nCL) {
        Y[(i-1)*nCL+j] <- rnorm(1,c(1,X[(i-1)*nCL+j]) %*% beta[k,],sigma)
        index[(i-1)*nCL+j] <- i
    }
}
B <- list()
for (i in 1:nBirder) {
    B[[i]] <- which(index == i)
}

# generate testing data
vX <- matrix(sqrt(runif(nBirder*nCL, min=0.001, max=100)),ncol=1)
vY <- array(0,nBirder*nCL)
index <- array(0,nBirder*nCL)
for (i in 1:nBirder) {
    k <- sample(1:K,1,prob=pi)
    
    for (j in 1:nCL) {
        vY[(i-1)*nCL+j] <- rnorm(1,c(1,vX[(i-1)*nCL+j]) %*% beta[k,],sigma)
        index[(i-1)*nCL+j] <- i
    }
}
vB <- list()
for (i in 1:nBirder) {
    vB[[i]] <- which(index == i)
}

# plot 
plot(X,Y)

# learning
lambda <- 0.01 # 0.01
model <- RandomRestartEM.SAC(K,X,Y,B,lambda,2) 
cat("Output true model:\n")
print(beta)
print(pi)
print(sigma)

# learn MSAC
Ks <- 1:4
models <- list()
results <- data.frame(K=array(0,length(Ks)),
                      tAvgLL=array(0,length(Ks)),
                      vAvgLL=array(0,length(Ks)))
bestAvgLL <- -Inf
besetK <- 1
for (K in Ks) {
    # train SACs
    model <- RandomRestartEM.SAC(K,X,Y,B,lambda,2) 
    models[[K]] <- model
    
    # compute log-likelihood
    AvgLL <- ComputeAvgLL.SAC(model$beta,model$pi,model$sigma,vX,vY,vB,lambda)
    cat("K is",K,"AvgLL is",AvgLL,"\n")
    if (bestAvgLL < AvgLL) {
        bestAvgLL <- AvgLL
        bestK <- K
    }
    
    results$K[K] <- K
    results$tAvgLL[K] <- model$AvgLL
    results$vAvgLL[K] <- AvgLL
} # K

cat("Best K is",bestK,"\n")
print(results)



