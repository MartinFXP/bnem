#' applies "inverse" absorption law to a disjuncitve normal form
#' @param bString a disjunctive normal form or binary vector according to model
#' @param model model for respective binary vector
#' @author Martin Pirkl
#' @return bString after "inverse" absorption law
#' @export
#' @examples
#' graph <- c("A+B=C", "A=C")
#' absorptionII(graph)
absorptionII <-
function(bString, model=NULL) {
  if (is.null(model)) {
    graph <- bString
    nodes <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(graph, "=")), "\\+"))))
  } else {
    graph <- model$reacID[which(bString == 1)]
    nodes <- model$namesSpecies
  }
  for (i in graph) {
    players <- unlist(strsplit(gsub("=.*", "", i), "\\+"))
    target <- gsub(".*=", "", i)
    others <- nodes[-which(nodes %in% c(players, target))]
    players2 <- gsub("!", "", players)
    change1 <- which(players == players2)
    change2 <- which(!(players == players2))
    if (length(change1) > 0) {
      others <- c(others, paste("!", players2[change1], sep = ""))
    }
    if (length(change2) > 0) {
      others <- c(others[-which(others %in% players2[change2])],  paste("\\+", players2[change2], sep = ""), paste("^", players2[change2], sep = ""))
    }
    targets <- intersect(grep(paste(paste("^", paste(players, collapse = "|^"), sep = ""), "|", paste("+", paste(players, collapse = "|+"), sep = ""), sep = ""), graph), grep(paste("=", target, sep = ""), graph))
    toomuch <- which(targets %in% grep(paste(others, collapse = "|"), graph))
    if (length(toomuch) > 0) {
      targets <- targets[-toomuch]
    }
    if (length(targets) > 1) {
      targets <- targets[-which(targets %in% which(graph %in% i))]
      if (is.null(model)) {
        if (sum(bString %in% graph[targets]) > 0) {
          bString <- bString[-which(bString %in% graph[targets])]
        }
      } else {
        bString[which(model$reacID %in% graph[targets])] <- 0
      }
    }
  }
  return(bString)
}
