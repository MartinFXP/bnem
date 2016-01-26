reduceGraph <- function(bString, model, CNOlist) {
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

## reduceGraph2 <- function(bString, model) {
##   graph <- model$reacID
##   for (i in graph) {
##     part <- unlist(strsplit(i, "="))
##     part <- c(unlist(strsplit(part[1], "\\+")), part[2])
##   }
##   ##if (length(colnames(CNOlist@stimuli)) == sum(stimuli %in% unique(part))) { # what is that for ?
##   ## handles all gatesizes
##   ## uses the fact that the gates in model$reacID are ordered by gatesize (really?)
##   parentsList <- list()
##   gates <- model$reacID
##   gatesList <- list()
##   if (length(gates) > 0) {
##     ## do pairwise reduction (kind of like an "pre"-extension to sperner): # that does not work, man!
##     ## for (i in which(bString == 1)) {
##     ##   for (j in which(bString == 1)) {
##     ##     if (i %in% j) { next() }
##     ##     edge1 <- unlist(strsplit(gates[i], "="))
##     ##     edge2 <- unlist(strsplit(gates[j], "="))
##     ##     if (!(edge1[2] %in% edge2[2])) { next() }
##     ##     output <- edge1[2]
##     ##     edge11 <- unlist(strsplit(edge1[1], "\\+"))
##     ##     edge22 <- unlist(strsplit(edge2[1], "\\+"))
##     ##     inputs <- intersect(edge11, edge22)
##     ##     if (length(inputs) == 0) { next() }
##     ##     bString[which(gates %in% paste(paste(inputs, collapse = "+"), output, sep = "="))] <- 1
##     ##   }
##     ## }
##     ## sperner reduction
##     for (i in 1:length(gates)) {
##       family <- unlist(strsplit(gates[i], "="))
##       child <- family[2]
##       parents <- unlist(strsplit(family[1], "\\+"))
##       gatesList[[i]] <- list()
##       gatesList[[i]]$parents <- parents
##       gatesList[[i]]$child <- child
##     }
##     for (i in 1:length(gates)) {
##       if (bString[i] == 1) {
##         for (j in 1:length(gatesList)) {
##           if (i == j) { next() }
##           if (sum(gatesList[[i]]$parents %in% gatesList[[j]]$parents == T) == length(gatesList[[i]]$parents) & gatesList[[i]]$child == gatesList[[j]]$child) {
##             if (bString[j] == 1) {
##               ##print(paste("removing: ", model$reacID[j], sep = ""))
##               bString[j] <- 0
##             }
##           }
##         }
##       }
##     }
##     ## delete non-signalling edges: does that even work? seems so, but lets confine it to a "small" number of gates
##     ## if (sum(bString == 1) < 50) { # too many?
##     ##   topNodes <- "start"
##     ##   while(!is.null(topNodes)) {
##     ##     usedGates <- gates[which(bString == 1)]
##     ##     topNodes <- NULL
##     ##     for (i in model$namesSpecies) {
##     ##       if (length(grep(paste(".*=", i, sep = ""), usedGates)) > 0 | i %in% colnames(CNOlist@stimuli)) {
##     ##         next()
##     ##       } else {
##     ##         gateTypeI <- grep(paste(".*\\+", i, ".*=", sep = ""), gates)
##     ##         if (length(gateTypeI) > 0) {
##     ##           if (1 %in% bString[gateTypeI]) {
##     ##                                     #print(paste("removing: ", gates[gateTypeI], sep = ""))
##     ##             bString[gateTypeI] <- 0
##     ##           }
##     ##         }
##     ##         gateTypeII <- grep(paste("^", i, ".*=", sep = ""), gates)
##     ##         if (length(gateTypeII) > 0) {
##     ##           if (1 %in% bString[gateTypeII]) {
##     ##                                     #print(paste("removing: ", gates[gateTypeII], sep = ""))
##     ##             bString[gateTypeII] <- 0
##     ##           }
##     ##         }
##     ##         if (length(intersect(c(grep(paste(".*\\+", i, ".*=", sep = ""), gates), grep(paste("^", i, ".*=", sep = ""), gates)), which(bString == 1))) > 0) {
##     ##           topNodes <- c(topNodes, i)
##     ##         }
##     ##       }
##     ##     }
##     ##   }
##   }
##   ## based on cnolist and the combination of perturbations one could even further reduce by gates that are never used for signaling (e.g. an and gate with input from two stimulations that are never combined in an experiment is not used for signaling (or at least can not be validated) and can be removed)
##   return(bString)
## }
