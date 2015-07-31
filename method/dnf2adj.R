dnf2adj <- function(dnf) {
  if (length(dnf) == 0) {
    return(NULL)
  } else {
    nodes <- character()
    for (i in dnf) {
      tmp <- unlist(strsplit(i, "="))
      nodes <- c(nodes, tmp[2], unlist(strsplit(tmp[1], "\\+")))
    }
    nodes <- unique(gsub("!", "", nodes))
    adjmat <- matrix(0, length(nodes), length(nodes))
    colnames(adjmat) <- nodes
    rownames(adjmat) <- nodes
    for (i in gsub("!", "", dnf)) {
      tmp <- unlist(strsplit(i, "="))
      child <- tmp[2]
      parents <- unlist(strsplit(tmp[1], "\\+"))
      adjmat[which(rownames(adjmat) %in% parents), which(colnames(adjmat) %in% child)] <- 1
    }
    diag(adjmat) <- 1
    stop <- FALSE
    cons <- c(TRUE, rep(FALSE, (length(adjmat) - 1)))
    while(!stop) {
      adjmat <- adjmat%*%adjmat
      if (all(cons == (adjmat != 0))) {
        stop <- TRUE
      } else {
        cons <- (adjmat != 0)
      }
    }
    adjmat[adjmat > 1] <- 1
    return(adjmat)
  }
}
