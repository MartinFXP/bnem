simulateStatesRecursiveIII <- function(CNOlist, model, bString, NEMlist = NULL) {
  require(matrixStats)
  getState <- function(CNOlist, node, signalStates, graph, children = NULL, NEMlist = NULL, bString, model) {
    graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
    if (length(graphCut) == 0) {
      signalStates[, node] <- 0
    } else {
      sop <- numeric(nrow(signalStates))
      children2 <- gsub("!", "", children)
      for (i in graphCut) {
        parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
        pob <- rep(1, nrow(signalStates))
        for (j in parents) {
          j2 <- gsub("!", "", j)
          if (sum(is.na(signalStates[, j2]) == T) == length(signalStates[, j2])) {
            if (j2 %in% children2) {
              print(length(bString))
              bString[which(model$reacID %in% i)] <- 0
              print(length(bString))
              graphCut <- graphCut[-which(graphCut %in% i)]
            } else {
              if (j %in% j2) {
                node2 <- node
                add1 <- 0
              } else {
                node2 <- paste("!", node, sep = "")
                add1 <- 1
              }
              get.state <- getState(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)), NEMlist, bString, model)
              signalStates <- get.state$signalStates
              bString <- get.state$bString
              graph <- get.state$graph
              if (add1 == 0) {
                pobMult <- signalStates[, j2]
              } else {
                pobMult <- add1 - signalStates[, j2]
              }
              #pobMult <- abs(add1 - signalStates[, j2])
              pobNA <- numeric(length(pob))
              pobNA[is.na(pob)] <- 1
              pobNA[is.na(pobMult)] <- 1
              pobMult[is.na(pobMult)] <- 1
              pob[is.na(pob)] <- 1

              pobMult[pobMult == -1] <- 0
              pob <- pob*pobMult
              
              pobNA[which(pob == 0)] <- 0
              pob[which(pobNA > 0)] <- NA
            }
          } else {
            if (j %in% j2) {
              add1 <- 0
            } else {
              add1 <- 1
            }
            if (add1 == 0) {
              pobMult <- signalStates[, j2]
            } else {
              pobMult <- add1 - signalStates[, j2]
            }
            pobNA <- numeric(length(pob))
            pobNA[is.na(pob)] <- 1
            pobNA[is.na(pobMult)] <- 1
            pobMult[is.na(pobMult)] <- 1
            pob[is.na(pob)] <- 1
            
            pobMult[pobMult == -1] <- 0
            pob <- pob*pobMult
            
            pobNA[which(pob == 0)] <- 0
            pob[which(pobNA > 0)] <- NA
          }
          if (max(pob, na.rm = T) == 0) { break() }
        }
        pobNA <- numeric(length(pob))
        pobNA[is.na(pob)] <- 1
        pobNA[is.na(sop)] <- 1
        pob[is.na(pob)] <- 0
        sop[is.na(sop)] <- 0
        sop <- sop + pob
        pobNA[which(sop > 0)] <- 0
        sop[which(pobNA > 0)] <- NA
        if (min(sop, na.rm = T) > 0) { break() }
      }
      sop[sop > 0] <- 1
      if (node %in% colnames(CNOlist@inhibitors)) {
        sop <- sop*(1 - CNOlist@inhibitors[, node]) # normal boolean way
      }
      if (node %in% colnames(CNOlist@stimuli)) {
        sop <- max(sop, CNOlist@stimuli[, node]) # normal boolean way
      }
      signalStates[, node] <- sop
    }
    return(list(signalStates = signalStates, bString = bString, graph = graph))
  }
  bString <- reduceGraph(bString, CNOlist, model)
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))]) 
  graph <- model$reacID[which(bString == 1)]
  stimuliStates <- CNOlist@stimuli
  if (!is.null(NEMlist$signalStates)) {
    signalStates <- NEMlist$signalStates
  } else {
    signalStates <- matrix(NA, nrow = nrow(CNOlist@signals[[2]]), ncol = length(inhibitors))
    rownames(signalStates) <- rownames(CNOlist@signals[[2]])
    colnames(signalStates) <- inhibitors
    signalStates <- cbind(stimuliStates, signalStates)
  }
  for (k in inhibitors) {
    if (sum(is.na(signalStates[, k]) == T) == length(signalStates[, k])) { ## this might lead to problems if there are signals that are neither inhibitors nor stimuli !
      get.state <- getState(CNOlist = CNOlist, node = k, signalStates = signalStates, graph = graph, children = NULL, NEMlist, bString, model)
      signalStates <- get.state$signalStates
      bString <- get.state$bString
    }
  }
  return(list(signalStates = signalStates, bString = bString))
}
