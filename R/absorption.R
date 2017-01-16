#' applies absorption law to a disjuncitve normal form
#' @param bString a disjunctive normal form or binary vector according to model
#' @param model model for respective binary vector
#' @author Martin Pirkl
#' @return bString after absorption law
#' @export
#' @examples
#' graph <- c("A+B=C", "A=C")
#' absorption(graph)
absorption <-
function(bString, model=NULL) {
  if (is.null(model)) {
    graph <- bString
  } else {
    graph <- model$reacID[which(bString == 1)]
  }
  for (i in graph) {
    targets <- grep(paste("(?=.*", gsub("\\+", ")(?=.*", gsub("=", ")(?=.*=", i)), ")", sep = ""), graph, perl = TRUE)
    toomuch <- grep(paste("!", gsub("\\+", "|!", gsub("=.*", "", i)), "", sep = ""), graph[targets])
    if (length(toomuch) > 0) {
      targets <- targets[-grep(paste("!", gsub("\\+", "|!", gsub("=.*", "", i)), "", sep = ""), graph[targets])]
    }
    if (length(targets) > 1) {
      targets <- targets[-which(targets == which(graph %in% i))]
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
