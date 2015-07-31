## dnf2adj <- function(dnf, nodes) {
##   adj <- matrix(0, length(nodes), length(nodes))
##   colnames(adj) <- nodes
##   rownames(adj) <- nodes
##   for (i in 1:length(nodes)) {
##     for (j in 1:length(nodes)) {
##       a <- nodes[i]
##       b <- nodes[j]
##       if (length(grep(paste(".*", a, ".*=", b, sep = ""), dnf)) > 0) {
##         adj[a, b] <- 1
##       }
##       if (length(grep(paste(".*", b, ".*=", a, sep = ""), dnf)) > 0) {
##         adj[b, a] <- 1
##       }
##     }
##   }
##   diag(adj) <- 1
##   return(adj)
## }

isDag <- function(graph = NULL, bString = 0, model = NULL) {
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

isDagOld <- function(bString, model) {
  dag <- TRUE
  dnf <- model$reacID[bString == 1]
  top <- model$namesSpecies[-which(model$namesSpecies %in% gsub(".*=", "", dnf))]
  layers <- list()
  layers[[1]] <- top
  count <- 1
  check <- ""
  while(length(top) > 0) {
    count <- count + 1
    top <- unique(gsub(".*=", "", dnf[grep(paste(paste(".*", layers[[(count-1)]], ".*=.*", sep = ""), collapse = "|"), dnf)]))
    if (sum(top %in% unlist(layers)) > 0) {
      layers[[(count-1)]] <- layers[[(count-1)]][-which(layers[[(count-1)]] %in% top)]
      #top <- top[-which(top %in% unlist(layers))]
    }
    if (length(top) > 0) {
      layers[[count]] <- sort(top)
    }
    if (all(unique(unlist(layers)) %in% check)) {
      break()
    }
    check <- unique(unlist(layers))
  }
  ranking <- unlist(layers)         
  adj <- dnf2adj(dnf)
  if (dim(adj)[1] > 1) {
    ranking <- ranking[which(ranking %in% rownames(adj))]
    adj <- adj[ranking, ranking]
    edges <- which(adj == 1, arr.ind = T)
    if (length(edges) != 0) {
      feedback <- NULL
      if (is.null(dim(edges))) {
        if (edges[1] > edges[2]) {
          feedback <- edges
        }
      } else {
        feedback <- edges[which(edges[, 1] > edges[, 2]), ]
      }
      if (length(feedback) > 1) {
        dag <- FALSE
      }
    }
  }
  return(dag)
}
