simulateStatesDAG <- function(CNOlist, model, bString, NEMlist = NULL) { ## does not work !!!
  signalStates <- matrix(0, nrow(CNOlist@cues), ncol(CNOlist@cues))
  colnames(signalStates) <- colnames(CNOlist@cues)
  rownames(signalStates) <- rownames(CNOlist@cues)
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- colnames(CNOlist@inhibitors)
  bString <- reduceGraph(bString, model)
  graph <- model$reacID[which(bString == 1)]
  hierarchy <- getHierarchy(graph)
  for (i in 1:length(hierarchy)) {
    for (j in 1:length(hierarchy[[i]])) {
      node <- hierarchy[[i]][j]
      if (node %in% stimuli) {
        signalStates[, which(colnames(signalStates) %in% node)] <- CNOlist@stimuli[, which(colnames(CNOlist@stimuli) %in% node)]
      } else {
        if (length(grep(paste("=", node, sep = ""), graph)) != 0) {
          sub.graph <- graph[grep(paste("=", node, sep = ""), graph)]
          sop <- numeric(nrow(CNOlist@cues))
          for (edge in sub.graph) {
            parents <- gsub(paste("=", node, sep = ""), "", unlist(strsplit(edge, "\\+")))
            pob <- rep(1, nrow(CNOlist@cues))
            for (parent in parents) {
              if (parent %in% gsub("!", "", parent)) {
                pob <- pob*signalStates[, which(colnames(signalStates) %in% parent)]
              } else {
                pob <- pob*(1-signalStates[, which(colnames(signalStates) %in% parent)])
              }
            }
            sop <- sop + pob
          }
          if (node %in% inhibitors) {
            sop <- sop*(1-CNOlist@inhibitors[, which(colnames(CNOlist@inhibitors) %in% node)])
          }
          sop[which(sop > 1)] <- 1
          signalStates[, which(colnames(signalStates) %in% node)] <- sop
        }
      }
    }
  }       
  return(signalStates)
}
