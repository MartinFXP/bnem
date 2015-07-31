makeGraphMat <- function(stimuli, inhibitors) {
  nodes <- c(stimuli, inhibitors)
  orGates <- matrix(0, length(nodes), length(inhibitors))
  colnames(orGates) <- inhibitors
  rownames(orGates) <- nodes
  andGates <- matrix(0, ((length(nodes)*(length(nodes) - 1))/2), length(inhibitors))
  colnames(andGates) <- inhibitors
  rownames(andGates) <- as.character(1:nrow(andGates))
  remember <- character()
  count <- 0
  for (i in nodes) {
    for (j in nodes) {
      if (paste(sort(c(i,j)), collapse = "+") %in% remember | i %in% j) { next() }
      remember <- c(remember, paste(sort(c(i,j)), collapse = "+"))
      count <- count + 1
      rownames(andGates)[count] <- paste(sort(c(i,j)), collapse = "+")
    }
  }
  return(list(orGates = orGates, andGates = andGates))
}
