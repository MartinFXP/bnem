makeSif <- function(stimuli, inhibitors, file = NULL) {
  if (is.null(file)) { stop("you have to specify a filename") }
  sifMatrix <- numeric()
  for (i in c(stimuli, inhibitors)) {
    for (j in inhibitors) {
      if (i %in% j) { next() }
      sifMatrix <- rbind(sifMatrix, c(i, 1, j))
      sifMatrix <- rbind(sifMatrix, c(i, -1, j))
    }
  }
  write.table(sifMatrix, file = file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
