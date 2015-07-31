simulateStatesRecursiveII <- function(CNOlist, model, bString) {
  require(matrixStats)
  bString <- reduceGraph(bString, CNOlist, model)
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))]) 
  graph <- model$reacID[which(bString == 1)]
  simulatedStates <- CNOlist@signals[[1]]
  simulatedStates[1:length(simulatedStates)] <- 0
  rownames(simulatedStates) <- rownames(CNOlist@signals[[1]])
  colnames(simulatedStates) <- colnames(CNOlist@signals[[1]])
  stimuliStates <- CNOlist@stimuli
  signalStates <- matrix(NA, nrow = nrow(CNOlist@inhibitors), ncol = length(inhibitors))
  rownames(signalStates) <- rownames(CNOlist@signals[[1]])
  colnames(signalStates) <- inhibitors
  signalStates <- cbind(stimuliStates, signalStates)
  for (k in inhibitors) {
    if (sum(is.na(signalStates[, k]) == T) > 0) {
      signalStates <- getStateII(CNOlist = CNOlist, node = k, signalStates = signalStates, graph = graph, children = NULL)
    }
  }
  #signalStates <- signalStates[, model$namesSpecies]
  signalStates[is.infinite(signalStates)] <- NA
  return(signalStates)
}

getStateII <- function(CNOlist, node, signalStates, graph, children = NULL) {
  graphCut <- graph[grep(paste("=", node, sep = ""), graph)]
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
        if (sum(is.na(signalStates[, j2]) == T) > 0) {
          if (j2 %in% children2) {
            if (j2 %in% j) {
              node2 <- node
              add1 <- 0
            } else {
              node2 <- paste("!", node, sep = "")
              add1 <- 1
            }
            if ((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2 != ceiling((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2)) {
              subGraph <- graph[-grep(paste(".*", node, ".*=.*", sep = ""), graph)]
              subResult <- getStateII(CNOlist = CNOlist, node = node, signalStates = signalStates, graph = subGraph, children = NULL)
              pobMult <- abs(add1 - subResult[, node])
              subResult <- signalStates
              subResult[, node] <- pobMult
              subGraph2 <- graph[-grep(paste(".*", j2, ".*=.*", sep = ""), graph)]
              subResult2 <- getStateII(CNOlist = CNOlist, node = j2, signalStates = subResult, graph = subGraph2, children = NULL)
              pobMult2 <- abs(add1 - subResult2[, j2])
              pobMult[which(pobMult != pobMult2)] <- -1
              if (node %in% colnames(CNOlist@inhibitors)) {
                pob <- pob*pobMult*(1 - CNOlist@inhibitors[, node])
              } else {
                pob <- pob*pobMult
              }
              pob[pob == -1] <- -Inf
            } else {
              subGraph <- graph[-grep(paste(".*", j2, ".*=.*", sep = ""), graph)]
              subResult <- getStateII(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = subGraph, children = NULL)
              pobMult <- abs(add1 - subResult[, j2])
              if (node %in% colnames(CNOlist@inhibitors)) {
                pob <- pob*pobMult*(1 - CNOlist@inhibitors[, node])
              } else {
                pob <- pob*pobMult
              }
            }
          } else {
            if (j %in% j2) {
              node2 <- node
              signalStates <- getStateII(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)))
              if (node %in% colnames(CNOlist@inhibitors)) {
                pob <- pob*signalStates[, j2]*(1 - CNOlist@inhibitors[, node])
              } else {
                pob <- pob*signalStates[, j2]
              }
            } else {
              node2 <- paste("!", node, sep = "")
              signalStates <- getStateII(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)))
              if (node %in% colnames(CNOlist@inhibitors)) {
                pob <- pob*(1 - signalStates[, j2])*(1 - CNOlist@inhibitors[, node])
              } else {
                pob <- pob*(1 - signalStates[, j2])
              }
            }
          }
        } else {
          if (j %in% j2) {
            if (node %in% colnames(CNOlist@inhibitors)) {
              pob <- pob*signalStates[, j2]*(1 - CNOlist@inhibitors[, node])
            } else {
              pob <- pob*signalStates[, j2]
            }
          } else {
            if (node %in% colnames(CNOlist@inhibitors)) {
              pob <- pob*(1 - signalStates[, j2])*(1 - CNOlist@inhibitors[, node])
            } else {
              pob <- pob*(1 - signalStates[, j2])
            }
          }
        }
        if (max(pob, na.rm = T) == 0) { break() }
      }
      sop <- sop + pob
      if (min(sop, na.rm = T) > 0) { break() }
    }
    sop[is.infinite(sop)] <- -Inf
    sop[sop > 0] <- 1
    signalStates[, node] <- sop
  }
  return(signalStates)
}
