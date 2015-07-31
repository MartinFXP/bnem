resolveTrue <- function(bString, new, model) {
  graph <- model$reacID[which(bString == 1)]
  new <- model$reacID[new]
  graph <- graph[-which(graph %in% new)]
  input <- unlist(strsplit(new, "="))
  output <- input[2]
  input <- unlist(strsplit(input[1], "\\+"))
  for (i in input) {
    if (i %in% gsub("!", "", i)) {
      if (length(grep(paste("!", i, ".*=", output, sep = ""), graph)) > 0) {
        graph <- graph[-grep(paste("!", i, ".*=", output, sep = ""), graph)]
      }
    } else {
      if (length(grep(paste("^", i, ".*=", output, "|\\+", i, ".*=", output, sep = ""), graph)) > 0) {
        graph <- graph[-grep(paste("^", i, ".*=", output, "|\\+", i, ".*=", output, sep = ""), graph)]
      }
    }
  }
  bString <- bString*0
  bString[which(model$reacID %in% c(new, graph))] <- 1
  return(bString)
}
