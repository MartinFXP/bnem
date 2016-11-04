complement <- function(bString, new, model) {
  graph <- model$reacID[which(bString == 1)]
  new2 <- gsub("!", "", new)
  if (new %in% new2) {
    delete <- grep(paste("!", new, sep = ""), graph)
  } else {
    delete <- grep(paste("^", new, "|\\+", new, sep = ""), graph)
  }
  graph <- graph[-delete]
  bString <- which(model$reacID %in% graph)
  return(bString)
}
    
