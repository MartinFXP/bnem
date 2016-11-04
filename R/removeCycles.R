removeCycles <-
function(bString, model, dnf = NULL) {
  if (is.null(dnf)) {
    if (any(bString != 0)) {
      graph <- model$reacID[which(bString == 1)]
      adjmat <- abs(dnf2adj(graph))
      get.order <- apply(adjmat, 2, sum)
      adjmat <- adjmat[order(get.order), order(get.order)]
      cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
      while(any(adjmat[lower.tri(adjmat)] == 1) & dim(cycles)[1] > 0) {
        for (i in 1:nrow(cycles)) {
          bad.cycles <- grep(paste(".*", rownames(adjmat)[cycles[i, 1]], ".*=", colnames(adjmat)[cycles[i, 2]], sep = ""), model$reacID)
          if (length(bad.cycles) > 0 & any(bString[bad.cycles] == 1)) {
            bString[bad.cycles] <- 0
            graph <- model$reacID[which(bString == 1)]
            break()
          }
        }
        adjmat <- abs(dnf2adj(graph))
        get.order <- apply(adjmat, 2, sum)
        adjmat <- adjmat[order(get.order), order(get.order)]
        cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
      }
      return(bString)
    } else {
      return(bString)
    }
  } else {
    graph <- dnf
    adjmat <- abs(dnf2adj(graph))
    get.order <- apply(adjmat, 2, sum)
    adjmat <- adjmat[order(get.order), order(get.order)]
    cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
    while(any(adjmat[lower.tri(adjmat)] == 1) & dim(cycles)[1] > 0) {
      for (i in 1:nrow(cycles)) {
        bad.cycles <- grep(paste(".*", rownames(adjmat)[cycles[i, 1]], ".*=", colnames(adjmat)[cycles[i, 2]], sep = ""), graph)
        if (length(bad.cycles) > 0) {
          graph <- graph[-bad.cycles]
          break()
        }
      }
      adjmat <- abs(dnf2adj(graph))
      get.order <- apply(adjmat, 2, sum)
      adjmat <- adjmat[order(get.order), order(get.order)]
      cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
    }
    return(graph)
  }
}
