BP <- function(CNOlist, NEMlist, model, sizeFac = 10^10, NAFac = 1, parameters = list(cutOffs = c(0.5,0.5,0), scoring = c(1,1,1)), method = "s", lambda = 1, beta = 1, limit = 10^-10, iter.max = 10) {
  P <- array(0.5, c(length(model$reacID), ncol(NEMlist$fc), 2))
  cavities <- expand.grid(1:length(model$reacID), 1:ncol(NEMlist$fc))
  P.dist <- Inf
  iter <- 0
  while(P.dist > limit & iter < iter.max) {
    iter <- iter + 1
    print(paste("Starting iteration ", iter, sep = ""))
    print(paste("Marginal Probability Distance: ", P.dist, sep = ""))
    P.old <- P
    pb <- txtProgressBar(style=1)
    end <- nrow(cavities)
    for (i in 1:end) {
      setTxtProgressBar(pb, i/end)
      cavity <- cavities[i, ]
      p.nu <- P[-cavity[1, 1],-cavity[1, 2],]
      P.mu <- lambda*exp(-cbind(rep(0, dim(p.nu)[1]), rep(1, dim(p.nu)[1])))*apply(p.nu, c(1,3), prod)
      P.mu <- P.mu/apply(P.mu, 1, sum)
      rand.string <- rep(0, length(model$reacID)-1)
      for (j in 1:(length(rand.string)-1)) {
        rand.string[j] <- sample(c(0,1), 1, prob = P.mu[j,])
      }
      #rand.string <- rbinom(length(model$reacID) - 1, 1, P.mu)
      if (cavity[1, 1] != 1 & cavity[1, 1] != length(model$reacID)) {
        test.string <- c(rand.string[1:(cavity[1, 1]-1)], 0, rand.string[cavity[1, 1]:length(rand.string)])
      }
      if (cavity[1, 1] == 1) {
        test.string <- c(0, rand.string[cavity[1, 1]:length(rand.string)])
      }
      if (cavity[1, 1] == length(model$reacID)) {
        test.string <- c(rand.string[1:(cavity[1, 1]-1)], 0)
      }
      ## do F on the subset mu (does not work with my scores then!!!):
      ## CNOlistTmp <- CNOlist
      ## NEMlistTmp <- NEMlist
      F <- computeScoreNemT1(CNOlist, model = NCNOcutCompExp, test.string, simList = NULL, indexList = NULL, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
      tmp.prod <- 1
      for (j in 1:nrow(P.mu)) {
        tmp.prod <- tmp.prod*P.mu[j, rand.string[j]+1]
      }
      #tmp.prod <- prod(P.mu[cbind(1:nrow(P.mu), rand.string+1)])
      P[cavity[1, 1],cavity[1, 2],1] <- beta*exp(-F)*tmp.prod
      if (cavity[1, 1] != 1 & cavity[1, 1] != length(model$reacID)) {
        test.string <- c(rand.string[1:(cavity[1, 1]-1)], 1, rand.string[cavity[1, 1]:length(rand.string)])
      }
      if (cavity[1, 1] == 1) {
        test.string <- c(1, rand.string[cavity[1, 1]:length(rand.string)])
      }
      if (cavity[1, 1] == length(model$reacID)) {
        test.string <- c(rand.string[1:(cavity[1, 1]-1)], 1)
      }
      F <- computeScoreNemT1(CNOlist, model = NCNOcutCompExp, test.string, simList = NULL, indexList = NULL, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
      tmp.prod <- 1
      for (j in 1:nrow(P.mu)) {
        tmp.prod <- tmp.prod*P.mu[j, rand.string[j]+1]
      }
      #tmp.prod <- prod(P.mu[cbind(1:nrow(P.mu), rand.string+1)])
      P[cavity[1, 1],cavity[1, 2],2] <- beta*exp(-F)*tmp.prod
      P[cavity[1, 1],cavity[1, 2],] <- P[cavity[1, 1],cavity[1, 2],]/sum(P[cavity[1, 1],cavity[1, 2],])
    }
    close(pb)
    P.dist <- sum((P-P.old)^2)
  }
  P <- exp(-cbind(rep(0, dim(P)[1]), rep(1, dim(P)[1])))*apply(P, c(1,3), prod)
  P <- P/apply(P, 1, sum)
  return(P)
}

## source("Boutros10.svn/method/scripts/cnopt.mod.R")

## system.time(marginals <- BP(CNOlist, NEMlist, model, sizeFac, parameters = parameters, method = method, iter.max = 1))

## m.prob <- cbind(marginals, TN)

## rownames(m.prob) <- model$reacID

## heatmapOP(m.prob, Colv = F, Rowv = T, breaks = seq(0,1,0.1))

## tmp <- numeric(10)

## for (i in 1:10000) {
  
## tmp <- tmp + rbinom(10, 1, seq(0.11,0.99, 0.09))

## }

## cbind(tmp/10000, seq(0.11,0.99, 0.09))

## test <- cbind(runif(10), runif(10, 1, 2))

## coords <- cbind(1:10, rep(c(1,2), 5))

## length(unlist(strsplit("==================================================================", "")))

MC <- function(CNOlist, NEMlist, model, sizeFac = 10^10, NAFac = 1, parameters = list(cutOffs = c(0.5,0.5,0), scoring = c(1,1,1)), method = "s", iter = 1000, parallel = NULL, parallel2 = NULL, verbose = FALSE) {
  
  B <- numeric(length(model$reacID))
   
  for (i in 1:iter) {
    
    bString <- sample(c(0,1), length(model$reacID), replace = T, prob = c(0.5,0.5))

    localString <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist, model=model, parameters=parameters, verbose = verbose, parallel=parallelList, parallel2=1, seeds=1, initSeed = bString, sizeFac = sizeFac, method = method)

  }

  return(B)

}

## source("Boutros10.svn/method/scripts/cnopt.mod.R"); system.time(marginals <- MC(CNOlist, NEMlist, model, sizeFac, parameters = parameters, method = method, iter = 10)); m.prob <- cbind(1-marginals, marginals, TN); rownames(m.prob) <- model$reacID; heatmapOP(m.prob, Colv = F, Rowv = T, breaks = seq(0,1,0.01)); m.prob[order(m.prob[, 1]), ]


## source("Boutros10.svn/method/scripts/cnopt.mod.R"); system.time(B <- MC(CNOlist, NEMlist, model, sizeFac, parameters = parameters, method = method, iter = 100)); B <- apply(B, 2, sum)/ncol(B); B <- B/sum(B); cbind(B[order(B)], TN[order(B)], model$reacID[order(B)])

