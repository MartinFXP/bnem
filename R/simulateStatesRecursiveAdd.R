simulateStatesRecursiveAdd <-
function(CNOlist, model, bString, NEMlist = NULL) { 
  require(matrixStats)
  require(Biobase)
  getStateAdd <- function(CNOlist, node, signalStates, graph, children = NULL, NEMlist = NULL) {
    graphCut <- graph[grep(paste("=", node, sep = ""), graph)]
    if (length(graphCut) == 0) {
      if (node %in% colnames(CNOlist@inhibitors)) {
        signalStates[, node] <- 0 - CNOlist@inhibitors[, node]
      } else {
        signalStates[, node] <- 0
      }
    } else {
      sop <- numeric(nrow(signalStates)) - 1
      children2 <- gsub("!", "", children)
      for (i in graphCut) {
        parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
        pob <- rep(1, nrow(signalStates))
        for (j in parents) {
          j2 <- gsub("!", "", j)
          if (sum(is.na(signalStates[, j2]) == T) == length(signalStates[, j2])) {
            if (j2 %in% children2) {
              if (j2 %in% j) {
                node2 <- node
                add1 <- 0
              } else {
                node2 <- paste("!", node, sep = "")
                add1 <- 1
              }
              if ((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2 != ceiling((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2)) {
                subGraph <- graph[-grep(paste(".*", node, ".*=.*", children2[length(children2)], sep = ""), graph)]
                subResult <- getStateAdd(CNOlist = CNOlist, node = node, signalStates = signalStates, graph = subGraph, children = NULL, NEMlist)
                pobMult <- subResult[, node]
                subResult <- signalStates
                subResult[, node] <- pobMult
                subGraph2 <- graph[-which(graph %in% i)]
                subResult2 <- getStateAdd(CNOlist = CNOlist, node = j2, signalStates = subResult, graph = subGraph2, children = NULL, NEMlist)
                if (add1 == 0) {
                  pobMult2 <- subResult2[, j2]
                } else {
                  pobMult2 <- add1 - subResult2[, j2]
                }
                #pobMult2 <- (add1 - subResult2[, j2])*(-1)^(add1+1)
                pobMult2[pobMult2 == 2] <- 1
                
                pobNA <- numeric(length(pob))
                pobNA[is.na(pob)] <- 1
                pobNA[is.na(pobMult)] <- 1
                pobNA[which(pobMult != pobMult2)] <- 1
                pobMult[is.na(pobMult)] <- 1
                pobMult[which(pobMult != pobMult2)] <- 1
                pob[is.na(pob)] <- 1
                
                ##pobMult[pobMult == -1] <- 0
                pob <- rowMin(cbind(pob,pobMult))
                
                pobNA[which(pob == 0)] <- 0
                pob[which(pobNA > 0)] <- NA
              } else {
                signalStatesTemp <- signalStates
                if (!(j %in% j2)) {
                  if (j2 %in% colnames(CNOlist@inhibitors)) {
                    signalStatesTemp[, node] <- 1 - CNOlist@inhibitors[, j2] # how to change this if the 2+ level states model is used?
                  } else {
                    signalStatesTemp[, node] <- 1
                  }
                } else {
                  if (j2 %in% colnames(CNOlist@inhibitors)) {
                    signalStatesTemp[, node] <- 0 - CNOlist@inhibitors[, j2] # how to change this if the 2+ level states model is used?
                  } else {
                    signalStatesTemp[, node] <- 0
                  }
                }
                subGraph <- graph[-which(graph %in% i)]
                subResult <- getStateAdd(CNOlist = CNOlist, node = j2, signalStates = signalStatesTemp, graph = subGraph, children = NULL, NEMlist)
                
                if (add1 == 0) {
                  pobMult <- subResult[, j2]
                } else {
                  pobMult <- add1 - subResult[, j2]
                }
                #pobMult <- (add1 - subResult[, j2])*(-1)^(add1+1)
                pobMult[pobMult == 2] <- 1
                
                pobNA <- numeric(length(pob))
                pobNA[is.na(pob)] <- 1
                pobNA[is.na(pobMult)] <- 1
                pobMult[is.na(pobMult)] <- 1
                pob[is.na(pob)] <- 1
                
                ##pobMult[pobMult == -1] <- 0
                pob <- rowMin(cbind(pob,pobMult))
                
                pobNA[which(pob == 0)] <- 0
                pob[which(pobNA > 0)] <- NA
              }
            } else {
              if (j %in% j2) {
                node2 <- node
                add1 <- 0
              } else {
                node2 <- paste("!", node, sep = "")
                add1 <- 1
              }
              signalStates <- getStateAdd(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)), NEMlist)
              
              if (add1 == 0) {
                  pobMult <- signalStates[, j2]
                } else {
                  pobMult <- add1 - signalStates[, j2]
                }
                #pobMult <- (add1 - signalStates[, j2])*(-1)^(add1+1)
              pobMult[pobMult == 2] <- 1
              
              pobNA <- numeric(length(pob))
              pobNA[is.na(pob)] <- 1
              pobNA[is.na(pobMult)] <- 1
              pobMult[is.na(pobMult)] <- 1
              pob[is.na(pob)] <- 1

              ##pobMult[pobMult == -1] <- 0
              pob <- rowMin(cbind(pob,pobMult))
              
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
                #pobMult <- (add1 - signalStates[, j2])*(-1)^(add1+1)
            pobMult[pobMult == 2] <- 1
            
            pobNA <- numeric(length(pob))
            pobNA[is.na(pob)] <- 1
            pobNA[is.na(pobMult)] <- 1
            pobMult[is.na(pobMult)] <- 1
            pob[is.na(pob)] <- 1
            
            ##pobMult[pobMult == -1] <- 0
            pob <- rowMin(cbind(pob,pobMult))
            
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
        sop <- rowMax(cbind(sop,pob))
        pobNA[which(sop > 0)] <- 0
        sop[which(pobNA > 0)] <- NA
        if (min(sop, na.rm = T) > 0) { break() }
      }
      ##sop[sop > 0] <- 1
      if (node %in% colnames(CNOlist@inhibitors)) {
        ##sop <- sop*(1 - CNOlist@inhibitors[, node]) # normal boolean way
        sop <- sop - CNOlist@inhibitors[, node] # kind of additive way to get more than two levels of state
      }
      if (node %in% colnames(CNOlist@stimuli)) {
        ##sop <- max(sop, CNOlist@stimuli[, node]) # normal boolean way
        sop <- sop + CNOlist@stimuli[, node] # kind of additive way to get more than two levels of state
      }
      signalStates[, node] <- sop
    }
    return(signalStates)
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
      signalStates <- getStateAdd(CNOlist = CNOlist, node = k, signalStates = signalStates, graph = graph, children = NULL, NEMlist)
    }
  }
  #inhibitorStates <- signalStates[, inhibitors]
  #inhibitorStates[which((1 - CNOlist@inhibitors) == 0)] <- 0
  #signalStates[, inhibitors] <- inhibitorStates
  return(signalStates)
}
