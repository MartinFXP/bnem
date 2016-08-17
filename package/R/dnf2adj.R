dnf2adj <-
function(dnf, closed = FALSE) {
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
    for (i in dnf) {
      tmp <- unlist(strsplit(i, "="))
      child <- tmp[2]
      parents <- unlist(strsplit(tmp[1], "\\+"))
      for (j in parents) {
        if (gsub("!", "", j) %in% j) {
          adjmat[which(rownames(adjmat) %in% j), which(colnames(adjmat) %in% child)] <- 1
        } else {
          adjmat[which(rownames(adjmat) %in% gsub("!", "", j)), which(colnames(adjmat) %in% child)] <- -1
        }
      }
    }
    diag(adjmat) <- 1
    if (closed) {
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
      adjmat[adjmat < -1] <- -1
    }
    return(adjmat)
  }
}
