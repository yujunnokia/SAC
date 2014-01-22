# The finite mixture of Species Accumulation Curves model
#
# Author: Jun Yu
# Version: Oct 2013
#
# Notations:
#   beta: coefficients 
#   pi: prior probabilities over K groups
#   sigma: sigma^2 is the shared variance
#   X: N by F feature matrix where N is the number of data instances and F is the number of features
#   Y: N by 1 label matrix where N is the number of data instances
#   B: a list of birders' observation index in X and Y
#   lambda: regularization penalty
###############################################################################

# magic numbers
EPSILON <- 1e-10
PI <- 3.1415926
MAX_ITERATIONS <- 50
SQRT2PI <- sqrt(2*PI)

#
# Compute average log-likelihood
#
ComputeAvgLL.SAC <- function(beta,pi,sigma,X,Y,B,lambda) 
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    K        <- length(pi)
    
    AvgLL <- 0
    for (i in 1:nBirder) {
        Yi <- Y[B[[i]]]
        Xi <- cbind(1,X[B[[i]],])
        
        iAvgLL <- 0
        for (k in 1:K) {
            ikAvgLL <- pi[k] * mean( exp(-(c(Xi %*% beta[k,]) - Yi)^2/(2*sigma^2)) / (SQRT2PI*sigma) )
            iAvgLL <- iAvgLL + ikAvgLL
        } # k
        
        if (iAvgLL == Inf || iAvgLL == -Inf) { stop("iAvgLL is Inf.\n") }
        if (is.na(iAvgLL)) { stop("iAvgLL is NA.\n") }
        if (is.nan(iAvgLL)) { stop("iAvgLL is NaN.\n") }
        
        if (iAvgLL == 0) {
            AvgLL <- AvgLL - 1/EPSILON
        } else {
            AvgLL <- AvgLL + log(iAvgLL)
        }
    } # i
    
    # regularization
    for (k in 1:K) {
        AvgLL <- AvgLL - nBirder*0.5*lambda*sum(beta[k,2:ncol(beta)]^2)
    } # k
    
    if (is.na(AvgLL)) { stop("AvgLL is na...\n") }
    return(as.numeric(AvgLL))
}

#
# Compute log-likelihood
#
ComputeLL.SAC <- function(beta,pi,sigma,X,Y,B,lambda) 
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    K        <- length(pi)
    
    LL <- 0
    for (i in 1:nBirder) {
        Yi <- Y[B[[i]]]
        Xi <- cbind(1,X[B[[i]],])
        nCL <- length(Yi)
        
        iLL <- 0
        for (k in 1:K) {
            ikLL <- exp( log(pi[k]) - sum( (c(Xi %*% beta[k,]) - Yi)^2 )/(2*sigma^2) - nCL*(log(SQRT2PI)+log(sigma)))
            iLL <- iLL + ikLL
        } # k
        
        if (iLL == Inf || iLL == -Inf) { stop("iLL is Inf.\n") }
        if (is.na(iLL)) { stop("iLL is NA.\n") }
        if (is.nan(iLL)) { stop("iLL is NaN.\n") }
        
        if (iLL == 0) {
            LL <- LL - 1/EPSILON
        } else {
            LL <- LL + log(iLL)
        }
    } # i
    
    # regularization
    for (k in 1:K) {
        LL <- LL - nBirder*0.5*lambda*sum(beta[k,2:ncol(beta)]^2)
    } # k
    
    if (is.na(LL)) { stop("LL is na...\n") }
    return(as.numeric(LL))
}


#
# Compute expected joint log-likelihood
#
ComputeEJLL.SAC <- function(beta,pi,sigma,X,Y,B,Z,lambda) 
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    K  <- length(pi)
    
    LOGSQRT2PI <- log(sqrt(2*PI) * sigma)
    
    EJLL <- 0
    for (i in 1:nBirder) {
        Yi <- Y[B[[i]]]
        Xi <- cbind(1,X[B[[i]],])
        nCL <- length(Yi)
        
        iEJLL <- 0
        for (k in 1:K) {
            ikEJLL <- Z[i,k] * (log(pi[k]) - sum( (c(Xi %*% beta[k,]) - Yi)^2 )/(2*sigma^2) - nCL*(log(SQRT2PI)+log(sigma)))
            iEJLL <- iEJLL + ikEJLL
        }
        
        if (iEJLL == Inf || iEJLL == -Inf) { stop("iEJLL is Inf.\n") }
        if (is.na(iEJLL)) { stop("iEJLL is NA.\n") }
        if (is.nan(iEJLL)) { stop("iEJLL is NaN.\n") }
        
        EJLL <- EJLL + iEJLL
    } # i
    
    # regularization
    for (k in 1:K) {
        EJLL <- EJLL - nBirder*0.5*lambda*sum(beta[k,2:ncol(beta)]^2) 
    } # k
    
    return(as.numeric(EJLL))
}

#
# Compute expected membership of birders in the E-step
#
ComputeZ.SAC <- function(beta,pi,sigma,X,Y,B)
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    K        <- length(pi)
    
    Z <- matrix(1,nrow=nBirder,ncol=K)
    for (i in 1:nBirder) {
        Yi  <- Y[B[[i]]]
        Xi <- cbind(1,X[B[[i]],])
        nCL <- length(Yi)
        
        denom <- 0
        for (k in 1:K) {
            Z[i,k] <- exp(log(pi[k]) - sum( (c(Xi %*% beta[k,]) - Yi)^2 )/(2 * sigma^2) )
            denom <- denom + Z[i,k]
            
            if (is.nan(Z[i,k])) { stop("Z(",i,",",k,") is nan.\n") }
        } # k
        
        # if denominator is zero, find the group with the largest probability and set it to be 1.
        if (denom == 0) { 
            maxProb <- -Inf
            largestK <- 1
            for (k in 1:K) {
                curProb <- log(pi[k]) - sum( (c(Xi %*% beta[k,]) - Yi)^2 )/(2 * sigma^2)
                if (maxProb < curProb) {
                    maxProb <- curProb
                    largestK <- k
                }
            } # k
            Z[i,largestK] <- 1
            denom <- 1
        }
        
        Z[i,] <- Z[i,] / denom
    } # i 
    
    Z[Z == 0] <- EPSILON
    Z[Z == 1] <- 1-EPSILON
    Z <- Z / rowSums(Z)
    
    return(Z)
}


#
# Expectation-Maximization
#
EM.SAC <- function(K,X,Y,B,lambda) 
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    
    # in the first E-step, randomly assign probabilities for Z
    Z <- matrix(runif(nBirder*K),nrow=nBirder)
    Z <- Z / rowSums(Z)
    
    beta <- matrix(0,nrow=K,ncol=(1+nFeature))
    newBeta <- beta
    models <- list()
    
    iteration <- 1
    diffParams <- 1.0e10
    while (diffParams > 1.0e-10 && iteration <= MAX_ITERATIONS) {
        # M-step        
        # update pi
        pi <- colMeans(Z)
        
        # update beta
        for (k in 1:K) {
            weights <- array(0,nData)
            for (i in 1:nBirder) {
                weights[B[[i]]] <- Z[i,k]
            }
            
            data <- data.frame(Y=Y,X=X)
            model <- glm( Y ~ X, data=data, weights=weights, family=gaussian())
            # model <- lm(Y ~ X, weights=weights)
            beta[k,] <- model$coefficients
            models[[k]] <- model
        }
        
        # update sigma
        term <- 0
        for (i in 1:nBirder) {
            Yi  <- Y[B[[i]]]
            Xi <- X[B[[i]],]
            nCL <- length(Yi)
            
            for (k in 1:K) {
                term <- term + Z[i,k] * sum( (c(cbind(1,Xi) %*% beta[k,]) - Yi)^2 )
            } # k
        } # i
        sigma <- sqrt(term / nData)
        
        # E-step
        Z <- ComputeZ.SAC(beta,pi,sigma,X,Y,B)
        
        # compute the change of parameters in two consecutive iterations
        oldBeta <- newBeta
        newBeta <- beta
        diffParams <- sum((newBeta-oldBeta)^2) / length(oldBeta)
        
        # output trace
        if (iteration %% 5 == 0) { 
            newLL <- ComputeLL.SAC(beta,pi,sigma,X,Y,B,lambda)
            newAvgLL <- ComputeAvgLL.SAC(beta,pi,sigma,X,Y,B,lambda)
            newEJLL <- ComputeEJLL.SAC(beta,pi,sigma,X,Y,B,Z,lambda)
            cat("EM iteration:",iteration, 
                    "LL:",newLL,
                    "AvgLL:",newAvgLL, 
                    "EJLL:",newEJLL, 
                    "params change:",diffParams,"\n") 
        }
        
        iteration <- iteration + 1
    }
    cat("EM converges in",iteration,"iterations...\n")
    newLL <- ComputeLL.SAC(beta,pi,sigma,X,Y,B,lambda)
    newAvgLL <- ComputeAvgLL.SAC(beta,pi,sigma,X,Y,B,lambda)
    newEJLL <- ComputeEJLL.SAC(beta,pi,sigma,X,Y,B,Z,lambda)
    
    return(list(beta=beta,pi=pi,sigma=sigma,Z=Z,models=models,LL=newLL,AvgLL=newAvgLL,EJLL=newEJLL))
}

#
# Random restart for EM
#
RandomRestartEM.SAC <- function(K,X,Y,B,lambda,nRandomRestarts=2) 
{
    nData    <- nrow(X)
    nFeature <- ncol(X) 
    nBirder  <- length(B)
    
    # random restart
    bestResult <- NULL
    bestLL <- bestAvgLL <- bestEJLL <- -Inf
    for (i in 1:nRandomRestarts) {
        cat("---------------------------\n")
        cat("Random restart EM run",i,"\n")
        
        done <- FALSE
        while (!done) 
        {
            r <- try(result <- EM.SAC(K,X,Y,B,lambda))
            done <- !inherits(r, "try-error")
        }
        
        cat("Final LL:", result$LL, "\n")
        cat("Final AvgLL:", result$AvgLL, "\n")
        cat("Final EJLL:", result$EJLL, "\n")
        
#        # print out model parameters
#        print(result$beta)
#        print(result$pi)
#        print(result$sigma)
        
        if (result$AvgLL > bestAvgLL) {
            bestResult <- result
            bestAvgLL <- result$AvgLL
            bestLL   <- result$LL
            bestEJLL <- result$EJLL
        }
    }
    
    # re-order the groups by the slope (beta_1)
    ordering <- order(bestResult$beta[,2],decreasing=TRUE)
    bestResult$pi <- bestResult$pi[ordering]
    bestResult$beta <- bestResult$beta[ordering,]
    bestResult$Z <- bestResult$Z[,ordering]
    bestResult$models <- bestResult$models[ordering]
    
    cat("***************************\n")
    cat("Best LL:", bestLL, "\n")
    cat("Best AvgLL:", bestAvgLL, "\n")
    cat("Best EJLL:", bestEJLL, "\n")
    
    if (K == 1) { 
        bestResult$beta <- matrix(bestResult$beta,nrow=1)
        bestResult$Z <- matrix(bestResult$Z,ncol=1)
    }
    cat("Output learned model:\n")
    print(bestResult$beta)
    print(bestResult$pi)
    print(bestResult$sigma)
    
    bestResult$lambda <- lambda
    return(bestResult)
}
