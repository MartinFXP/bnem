makeDesign2 <- function(x, stimuli, inhibitors, batches = NULL, runs = NULL) {
  design <- numeric()
  designNames <- character()
  design2 <- numeric()
  designNames2 <- character()
  for (i in c(stimuli, inhibitors, batches, runs)) {
    tmp <- numeric(ncol(x))
    tmp2 <- numeric(ncol(x))
    tmp[grep(paste("^", i, "|_", i, sep = ""), colnames(x))] <- 1
    tmp2[grep(paste("!", i, sep = ""), colnames(x))] <- 1
    if (sum(tmp) != 0) {
      design <- cbind(design, tmp)
      designNames <- c(designNames, i)
    }
    if (sum(tmp2) != 0) {
      design2 <- cbind(design2, tmp2)
      designNames2 <- c(designNames2, i)
    }
  }
  colnames(design) <- designNames
  colnames(design2) <- designNames2
  return(list(stimuli=design, inhibitors=design2))
}
