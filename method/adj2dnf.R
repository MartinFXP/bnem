adj2dnf <- function(A) {

  dnf <- NULL
  
  for (i in 1:ncol(A)) {
    for (j in 1:nrow(A)) {
      if (i %in% j) { next() }
      if (A[i, j] == 1) {
        dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
      }
    }
  }

  dnf <- unique(dnf)
  
  return(dnf)

}

## the following two function are not mine but were given to me by Benedict Anchang or Katharina Meyer (cannot remember whom):

## Adjacency marix to graph

## graphNEL-object from adj-matrix
adj2graph <- function(adj.matrix) {
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


# graph to adjacency matrix

graph2adj <- function(gR) {
    adj.matrix <- matrix(0,
                       length(nodes(gR)),
                       length(nodes(gR))
                       )
    rownames(adj.matrix) <- nodes(gR)
    colnames(adj.matrix) <- nodes(gR)
    for (i in 1:length(nodes(gR))) {
    adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
    }

    return(adj.matrix)
  }


