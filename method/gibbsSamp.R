### the following does NOT do the marginals! But estimates them?:

gibbsSamp <- function(CNOlist, NEMlist2, model, approach="fc", parameters= c(cutOffs = c(0.5,0.5,0), scoring = c(0.25,0.5,2)), verbose = FALSE, seeds = 1, iterations = 1000, burnin = 1000, parallel=NULL, sizeFac = 1, draw = FALSE) {
  ## do the burn in phase:
  bString <- numeric(length(model$reacID))
  marg.probs <- matrix(0, burnin, length(model$reacID))
  for (i in 1:burnin) {
    start <- Sys.time()
    for (j in 1:length(model$reacID)) {
      if (verbose) {
        print(paste("burnin step ", i, " edge ", j, sep = ""))
      }
      bString[j] <- 0
      tmp <- computeScoreNemT1(CNOlist, model, bString, sizeFac=sizeFac, approach = approach, NEMlist = NEMlist2, parameters=parameters, tellme = 0)
      bString[j] <- 1
      if (absorption(bString, model)[j] == 0) {
        bString <- absorptionII(bString, model)
      }
      tmp2 <- computeScoreNemT1(CNOlist, model, bString, sizeFac=sizeFac, approach = approach, NEMlist = NEMlist2, parameters=parameters, tellme = 0)
      marg.probs[i, j] <- tmp2/(tmp+tmp2)
      if (marg.probs[i, j] >= 0.5) {
        bString[j] <- 0
      }
      marg.probs[i, j] <- sample(c(0,1), 1, prob = c(1 - tmp2/(tmp+tmp2), tmp2/(tmp+tmp2)))
      if (verbose) {
        print(paste("network: ", toString(bString), sep = ""))
        print(paste("score: ", min(c(tmp, tmp2)), sep = ""))
      }
    }
    if (verbose) {
      if (i > 2) {
        plot(apply(marg.probs[1:i, ], 2, mean), xlab = "edges", ylab = "edge probabilities")
        abline(h = 0.5, col = "red")
      }
    }
    end <- Sys.time()
    print(end-start)
  }
  if (draw) {
    plotDnf(model$reacID[as.logical(bString)])
  }
  real.marg.probs <- matrix(0, iterations, length(model$reacID))
  for (i in 1:iterations) {
    start <- Sys.time()
    for (j in 1:length(model$reacID)) {
      if (verbose) { print(paste("iteration step ", i, " edge ", j, sep = "")) }
      bString[j] <- 0
      tmp <- computeScoreNemT1(CNOlist, model, bString, sizeFac=sizeFac, approach = approach, NEMlist = NEMlist2, parameters=parameters, tellme = 0)
      bString[j] <- 1
      if (absorption(bString, model)[j] == 0) {
        bString <- absorptionII(bString, model)
      }
      tmp2 <- computeScoreNemT1(CNOlist, model, bString, sizeFac=sizeFac, approach = approach, NEMlist = NEMlist2, parameters=parameters, tellme = 0)
      real.marg.probs[i, j] <- tmp2/(tmp+tmp2)
      if (real.marg.probs[i, j] >= 0.5) {
        bString[j] <- 0
      }
      real.marg.probs[i, j] <- sample(c(0,1), 1, prob = c(1 - tmp2/(tmp+tmp2), tmp2/(tmp+tmp2)))
      if (verbose) {
        print(paste("network: ", toString(bString), sep = ""))
        print(paste("score: ", min(c(tmp, tmp2)), sep = ""))
      }
    }
    if (verbose) {
      if (i > 2) {
        plot(apply(real.marg.probs[1:i, ], 2, mean), xlab = "edges", ylab = "edge probabilities")
        abline(h = 0.5, col = "red")
      }
    }
    end <- Sys.time()
    print(end-start)
  }
  return(list(bString=bString, burnin=marg.probs, real=real.marg.probs))
}
