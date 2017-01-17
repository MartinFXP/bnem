#' @noRd
isDag <-
function(graph = NULL, bString = 0, model = NULL) {
  if (any(bString != 0)) {
    graph <- model$reacID[which(bString == 1)]
  }
  if (!is.null(graph)) {
    adjmat <- dnf2adj(graph)
    ##.order <- apply(adjmat, 1, sum)
    get.order2 <- apply(adjmat, 2, sum)
    adjmat <- adjmat[order(get.order2, decreasing = F), order(get.order2, decreasing = F)]
    if (all(adjmat[lower.tri(adjmat)] == 0)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}
