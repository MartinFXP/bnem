plotBinary <- function(bString, model, ...) {
  #gates <- model$reacID[which(bString == 1)]
  #for (i in 1:length(grep("\\+", gates))) {
  #  graphNodes <- c(graphNodes, paste("AND", i, sep = ""))
  #}
  #nodes <- character()
  #for (gate in gates) {
  #  temp <- unlist(strsplit(gate, split="="))
  #  graphNodes <- c(graphNodes, temp[2])
  #  temp <- gsub("!", "", unlist(strsplit(temp[1], split="\\+")))
  #  graphNodes <- c(graphNodes, temp)
  #}
  #graphNodes <- unique(graphNodes)
  #adjGraph <- matrix(0, nrow = length(graphNodes), ncol = length(graphNodes))
  #rownames(adjGraph) <- graphNodes
  #colnames(adjGraph) <- graphNodes
  #andGate <- 1
  #for (gate in gates) {
  #  temp <- unlist(strsplit(gate, split="="))
  #  inNode <- temp[2]
  #  temp <- gsub("!", "", unlist(strsplit(temp[1], split="\\+")))
  #  if (length(temp) == 1) {
  #    outNode <- temp
  #    adjGraph[which(rownames(adjGraph) == outNode), which(colnames(adjGraph) == inNode)] <- 1
  ModelCut <- model
  ModelCut$interMat <- ModelCut$interMat[, as.logical(bString)]
  ModelCut$notMat <- ModelCut$notMat[, as.logical(bString)]
  ModelCut$reacID <- ModelCut$reacID[as.logical(bString)]
  plotModel(ModelCut, ...)
}
