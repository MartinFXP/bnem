#' @noRd
adj2graph <-
function(adj.matrix) {
  V   <- rownames(adj.matrix)
  edL <- vector("list", length=nrow(adj.matrix))
  names(edL) <- V
  for (i in 1:nrow(adj.matrix)) {
    edL[[i]] <- list(edges=which(!adj.matrix[i,]==0),
                     weights=adj.matrix[i,!adj.matrix[i,]==0])
  }
  gR <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
  return(gR)
}
