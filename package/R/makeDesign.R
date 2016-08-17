makeDesign <-
function(x, stimuli, inhibitors, batches = NULL, runs = NULL) {
  design <- numeric()
  designNames <- character()
  for (i in c(stimuli, inhibitors, batches, runs)) {
    tmp <- numeric(ncol(x))
    tmp[grep(i, colnames(x))] <- 1
    if (sum(tmp) != 0) {
      design <- cbind(design, tmp)
      designNames <- c(designNames, i)
    }
  }
  colnames(design) <- designNames
  return(design)
}
