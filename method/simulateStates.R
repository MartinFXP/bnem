simulateStates <- function(CNOlist, model, bString) {
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))]) # colnames(CNOlist@inhibitors)
  graph <- model$reacID[which(bString == 1)]
  simulatedStates <- CNOlist@signals[[1]]
  simulatedStates[1:length(simulatedStates)] <- 0
  rownames(simulatedStates) <- rownames(CNOlist@signals[[1]])
  colnames(simulatedStates) <- colnames(CNOlist@signals[[1]])
  stimuliStates <- CNOlist@stimuli
  signalStates <- cbind((1 - CNOlist@inhibitors), matrix(1, nrow = nrow(CNOlist@inhibitors), ncol = (length(inhibitors) - ncol(CNOlist@inhibitors))))
  signalStates[which(signalStates == 1)] <- NA
  signalStates <- signalStates
  rownames(signalStates) <- rownames(CNOlist@signals[[1]])
  colnames(signalStates) <- c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))])
  signalStates <- cbind(stimuliStates, signalStates)
  gates <- list()
  parents <- list()
  for (i in inhibitors) {
    gates[[i]] <- graph[grep(paste("=", i, sep = ""), graph)]
    if (length(gates) == 0) { next() }
    for (j in 1:length(gates[[i]])) {
      if (length(grep("\\+", gates[[i]][j])) >= 1) {
        tmp <- unlist(strsplit(gates[[i]][j], "\\+"))
        tmp2 <- unlist(strsplit(tmp[length(tmp)], "="))
        parents[[i]] <- unique(c(parents[[i]], tmp[-length(tmp)], tmp2[1]))
      } else {
        tmp <- unlist(strsplit(gates[[i]][j], "="))
        parents[[i]] <- unique(c(parents[[i]], tmp[1]))
      }
    }
  }
  ## try du built adjacency matrix to get a good order of inhibitors to check
  nodes <- sort(gsub("!", "", unique(unlist(strsplit(unlist(strsplit(model$reacID, "\\+")), "=")))))
  andGates <- grep("\\+", model$reacID)
  orGates <- 1:length(model$reacID)
  orGates <- orGates[-andGates]
  orGates <- model$reacID[orGates]
  adjMat <- matrix(0, length(nodes), length(nodes))
  colnames(adjMat) <- nodes
  rownames(adjMat) <- nodes
  for (i in nodes) {
    tmp <- gsub("!", "", gsub(paste("=", i, sep = ""), "", orGates[grep(paste("=", i, sep = ""), orGates)]))
    if (length(tmp) > 0) {
      signs <- numeric(length(tmp)) + 1
      signs[grep("!", gsub(paste("=", i, sep = ""), "", orGates[grep(paste("=", i, sep = ""), orGates)]))] <- -1
      adjMat[tmp, i] <- signs
    }
  }
  require(expm)
  for (i in 2:length(nodes)) {
    adjMat <- adjMat + adjMat%^%i
  }
  RS <- rowSums(abs(adjMat))
  RSord <- order(RS, decreasing = T)
  adjMat <- adjMat[RSord, RSord]
  adjMat2 <- adjMat
  adjMat2[which(adjMat2 > 0)] <- 1
  adjMat2[which(adjMat2 < 0)] <- -1
  #print(adjMat2)
  inhibitorsOrd <- intersect(colnames(adjMat), inhibitors)
  #inhibitorsOrd <- inhibitors
  ## get states
  parents2 <- list()
  for (j in inhibitorsOrd) {
    sop <- numeric(nrow(CNOlist@stimuli))
    for (k in grep(paste("=", j, sep = ""), graph)) {
      tmp <- unlist(strsplit(graph[k], "\\="))
      tmp2 <- unlist(strsplit(tmp[1], "\\+"))
      pob <- rep(1, nrow(CNOlist@stimuli))
      for (l in tmp2) {
        if (length(grep("!", l)) == 0) {
          if (j %in% colnames(CNOlist@inhibitors)) {
            pob <- pob*signalStates[, l]*(1 - CNOlist@inhibitors[, j])
          } else {
            pob <- pob*signalStates[, l]
          }
        } else {
          if (j %in% colnames(CNOlist@inhibitors)) {
            pob <- pob*NOT(signalStates[, gsub("!", "", l)])*(1 - CNOlist@inhibitors[, j])
          } else {
            pob <- pob*NOT(signalStates[, gsub("!", "", l)])
          }
        }
      }
      sop <- cbind(sop, pob)
    }
    if (is.null(dim(sop))) {
      signalStates[, j] <- sop
    } else {
      signalStates[, j] <- apply(sop, 1, max)
    }
  }
  #signalStates <- signalStates[, colnames(CNOlist@inhibitors)]
  #signalStates <- cbind(CNOlist@stimuli, signalStates)
  signalStates <- signalStates[, model$namesSpecies]
  ## for (i in 1:nrow(CNOlist@stimuli)) {
  ##   count <- 0
  ##   while(sum(is.na(signalStates[i, ])) > 0 & count <= length(inhibitors)) {
  ##     count <- count + 1
  ##     for (j in inhibitorsOrd) {
  ##       if (is.null(parents[[j]])) {
  ##         signalStates[i, j] <- 0
  ##       } else {
  ##         if (sum(is.na(signalStates[i, gsub("!", "", parents[[j]])]) == TRUE) > 0) { next() }
  ##         statesTmp <- numeric()
  ##         for (k in 1:length(gates[[j]])) {
  ##           #if (length(grep("\\+", gates[[j]][k])) < 1) { next() }
  ##           parents2[[k]] <- character()
  ##           tmp <- unlist(strsplit(gates[[j]][k], "\\+"))
  ##           tmp2 <- unlist(strsplit(tmp[length(tmp)], "="))
  ##           parents2[[k]] <- unique(c(parents2[[k]], tmp[-length(tmp)], tmp2[1]))
  ##           parStatesTmp <- numeric()
  ##           for (l in parents2[[k]]) {
  ##             if (length(grep("!", l)) == 1) {
  ##               parStatesTmp <- c(parStatesTmp, NOT(signalStates[i, gsub("!", "", l)]))
  ##             } else {
  ##               parStatesTmp <- c(parStatesTmp, OR(signalStates[i, l]))
  ##             }
  ##           }
  ##           statesTmp <- c(statesTmp, min(parStatesTmp))
  ##         }
  ##         signalStates[i, j] <- max(statesTmp)
  ##       }
  ##     }
  ##   }
    #print(count)
  #}
  return(signalStates)
}

NOT <- function(x) { return(1 - x) }
OR <- function(x) { return(x) }
