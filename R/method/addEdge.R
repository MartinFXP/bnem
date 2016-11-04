addEdge <- function(edges, CNOlist, model, n = 100, full = FALSE) {
  sifMatrix <- numeric()
  graph <- model$reacID[-grep("\\+", model$reacID)]
  for (i in c(graph, edges)) {
    tmp2 <- unlist(strsplit(i, "="))
    if (gsub("!", "", tmp2[1]) %in% tmp2[1]) {
      sifMatrix <- rbind(sifMatrix, c(tmp2[1], 1, tmp2[2]))
    } else {
      sifMatrix <- rbind(sifMatrix, c(gsub("!", "", tmp2[1]), -1, tmp2[2]))
    }
  }
  write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  PKN2 <- readSIF("temp.sif")
  unlink("temp.sif")
  model2 <- preprocessing(CNOlist, PKN2, maxInputsPerGate=n)
  if (!full) {
    index <- which(!(model2$reacID %in% model$reacID) & !(model2$reacID %in% edges))
    model2$reacID <- model2$reacID[-index]
    model2$interMat[, -index]
    model2$notMat[, -index]
  }
  return(model2)
}
