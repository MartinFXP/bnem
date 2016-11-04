addDumpnodes <- function(CNOlist, names = c("D1", "D2")) {
  for (i in 1:length(CNOlist@signals)) {
    CNOlist@signals[[i]] <- cbind(CNOlist@signals[[i]], A = 0, B = 0)
    colnames(CNOlist@signals[[i]])[(ncol(CNOlist@signals[[i]]) - length(names) + 1):ncol(CNOlist@signals[[i]])] <- names
  }
  return(CNOlist)
}
      
