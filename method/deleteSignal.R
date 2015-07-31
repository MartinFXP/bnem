deleteSignal <- function(s, CNOlist) {
  CNOlist2 <- CNOlist
  for (i in s) {
    if (i %in% colnames(CNOlist2@signals[[1]])) {
      CNOlist2@signals[[1]] <- CNOlist2@signals[[1]][, -which(colnames(CNOlist2@signals[[1]]) %in% i)]
      CNOlist2@signals[[2]] <- CNOlist2@signals[[2]][, -which(colnames(CNOlist2@signals[[2]]) %in% i)]
    }
    if (i %in% colnames(CNOlist2@stimuli)) {
      CNOlist2@cues <- CNOlist2@cues[, -which(colnames(CNOlist2@cues) %in% i)]
      if (!is.null(dim(CNOlist2@stimuli[, -which(colnames(CNOlist2@stimuli) %in% i)]))) {
        CNOlist2@stimuli <- CNOlist2@stimuli[, -which(colnames(CNOlist2@stimuli) %in% i)]
      } else {
        CNOlist2@stimuli <- as.matrix(CNOlist2@stimuli[, -which(colnames(CNOlist2@stimuli) %in% i)])
      }
    }
    if (i %in% colnames(CNOlist2@inhibitors)) {
      CNOlist2@cues <- CNOlist2@cues[, -which(colnames(CNOlist2@cues) %in% i)]
      if (!is.null(dim(CNOlist2@inhibitors[, -which(colnames(CNOlist2@inhibitors) %in% i)]))) {
        CNOlist2@inhibitors <- CNOlist2@inhibitors[, -which(colnames(CNOlist2@inhibitors) %in% i)]
      } else {
        CNOlist2@inhibitors <- as.matrix(CNOlist2@inhibitors[, -which(colnames(CNOlist2@inhibitors) %in% i)])
      }
    }
  }
  return(CNOlist2)
}
