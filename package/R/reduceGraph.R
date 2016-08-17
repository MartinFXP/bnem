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
