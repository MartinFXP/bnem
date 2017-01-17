#' reduces the size of a graph if ossible to an equivalent sub-graph
#' @param bString binary vector indicating the sub-graph given a model
#' @param model model object for the whole graph space
#' @param CNOlist CNOlist object
#' @author Martin Pirkl
#' @return equivalent sub-graph denoted by bString
#' @export
#' @examples
#' library(bnem)
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- readSIF("temp.sif")
#' unlink('temp.sif')
#' model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
#' bString <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)
reduceGraph <-
function(bString, model, CNOlist) {
  if (any(bString != 0)) {
    stimuli <- colnames(CNOlist@stimuli)
    graph <- model$reacID[which(bString == 1)]
    tmp <- unlist(strsplit(graph, "="))
    tmp <- unlist(strsplit(tmp, "\\+"))
    tmp <- unique(gsub("!", "", tmp))
    for (i in tmp) {
      if (!(i %in% stimuli) & length(grep(paste("=", i, sep = ""), graph)) == 0) {
        ## get <- grep(paste("\\+", i, "|^", i, sep = ""), graph) # this is not good
        get <- grep(paste("^", i, sep = ""), graph) # this is better; more conservative
        if (length(get) > 0) {
          graph <- graph[-get]
        }
      }
    }
    bString <- numeric(length(bString))
    bString[which(model$reacID %in% graph)] <- 1
  }
  bString <- absorption(bString, model)
  return(bString)
}


