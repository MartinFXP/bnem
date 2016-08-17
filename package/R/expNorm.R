expNorm <-
function(x,  stimuli = NULL, inhibitors = NULL, batches, runs, cutoff) {
  design <- makeDesign(x, stimuli, inhibitors, c(batches, runs))
  # for every gene find all the batches in all the runs that have a significant change and normalize (pam or simple)
  # normalize batches that do not have a significant change according to a level cutoff or to control level
  normedX <- x*0
  for (run in runs) {
    for (batch in batches) {
      targetRows <- intersect(which(design[, run] == 1), which(design[, batch] == 1))
      if (length(targetRows) == 0) { next() }
      for (i in 1:nrow(normedX)) {
        if (max(x[i, targetRows]) - min(x[i, targetRows]) >= cutoff) {
          normedX[i, targetRows] <- simpleNorm(x[i, targetRows])
        }
      }
    }
  }
  cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors), collapse = "|"), colnames(design))], 1, sum)
  grepCtrl <- which(cuesSum == 0)
  for (run in runs) {
    for (batch in batches) {
      targetRows <- intersect(which(design[, run] == 1), which(design[, batch] == 1))
      if (length(targetRows) == 0) { next() }
      for (i in 1:nrow(normedX)) {
        if (max(x[i, targetRows]) - min(x[i, targetRows]) < cutoff) {
          normedX[i, targetRows] <- median(normedX[i, grepCtrl])
        }
      }
    }
  }
  return(normedX)
}
