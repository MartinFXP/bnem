checkCNOlist <- function(CNOlist) {
  if (dim(CNOlist@stimuli)[2] == 1) {
    name <- colnames(CNOlist@stimuli)
    CNOlist@stimuli <- as.matrix(CNOlist@stimuli[, order(colnames(CNOlist@stimuli))])
    colnames(CNOlist@stimuli) <- sort(name)
  }else {
    CNOlist@stimuli <- CNOlist@stimuli[, order(colnames(CNOlist@stimuli))]
  }
  if (dim(CNOlist@inhibitors)[2] == 1) {
    name <- colnames(CNOlist@inhibitors)
    CNOlist@inhibitors <- as.matrix(CNOlist@inhibitors[, order(colnames(CNOlist@inhibitors))])
    colnames(CNOlist@inhibitors) <- sort(name)
  } else {
    CNOlist@inhibitors <- CNOlist@inhibitors[, order(colnames(CNOlist@inhibitors))]
  }
  CNOlist@cues <- cbind(CNOlist@stimuli, CNOlist@inhibitors)
  for (i in 1:length(CNOlist@signals)) {
    CNOlist@signals[[i]] <- CNOlist@signals[[i]][, order(colnames(CNOlist@signals[[i]]))]
  }
  CNOlist@variances <- list()
  return(CNOlist)
}
