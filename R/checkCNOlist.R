#' @noRd
checkCNOlist <-
function(CNOlist) {
  if (dim(CNOlist@stimuli)[2] > 1) {
    CNOlist@stimuli <- CNOlist@stimuli[, order(colnames(CNOlist@stimuli))]
  }
  if (dim(CNOlist@inhibitors)[2] > 1) {
    CNOlist@inhibitors <- CNOlist@inhibitors[, order(colnames(CNOlist@inhibitors))]
  }
  CNOlist@cues <- cbind(CNOlist@stimuli, CNOlist@inhibitors)
  if (ncol(CNOlist@signals[[1]]) > 1) {
    for (i in 1:length(CNOlist@signals)) {
      CNOlist@signals[[i]] <- CNOlist@signals[[i]][, order(colnames(CNOlist@signals[[i]]))]
    }
  }
  CNOlist@variances <- list()
  return(CNOlist)
}
