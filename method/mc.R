mc <- function(CNOlist, NEMlist, model, sizeFac = 10^10, NAFac = 1, parameters = list(cutOffs = c(0.5,0.5,0), scoring = c(1,1,1)), method = "s", size = 1000, stop = 1000, max.draws = 10000, prob = c(0.5,0.5), n = 10, verbose = TRUE) {

  top <- matrix(0, size, length(model$reacID))

  score <- numeric(size)

  count <- 0

  iter <- 0

  while(iter < stop & count < max.draws+size) {

    if (count >= size) {

      draw <- sample(c(0,1), length(model$reacID), replace = T, prob = prob)
      
      cost <- computeScoreNemT1(CNOlist, model = NCNOcutCompExp, draw, simList = NULL, indexList = NULL, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)

      if (cost > 0) { cost <- 0 }

      if (cost < max(score)) {

        keep <- 1

      } else {
      
        itsprob <- c(0.999, 0.001)*c(1 - cost/max(score), cost/max(score))

        print(itsprob)
        
        keep <- sample(c(0,1), 1, prob = itsprob)

      }

      if (keep == 1) {

        top[which(score == max(score)), ] <- draw

        score[which(score == max(score))] <- cost

        iter <- 0

      } else {

        iter <- iter + 1

      }
      
      count <- count + 1

    } else {
      
      count <- count + 1
      
      draw <- sample(c(0,1), length(model$reacID), replace = T, prob = prob)
      
      cost <- computeScoreNemT1(CNOlist, model = NCNOcutCompExp, draw, simList = NULL, indexList = NULL, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)

      if (cost > 0) { cost <- 0 }

      top[count, ] <- draw

      score[count] <- cost

    }

    if (verbose) {
      
      print(paste("Sample: ", count, " - Stall: ", iter, " - Best Score: ", min(score), sep = ""))
      
    }

  }

  probs <- colSums(top)/size

  bStrings <- matrix(0, n, length(model$reacID))
  
  for (i in 1:n) {
    
    for (j in 1:length(model$reacID)) {
      
      bStrings[i, j] <- sample(c(0,1), 1, prob = c(1 - probs[j], probs[j]))
      
    }
    
  }
  
  return(bStrings)

}
