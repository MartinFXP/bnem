getHierarchy <- function(graph) {
  adj <- dnf2adj(graph)
  require(nem)
  adj2 <- transitive.reduction(abs(adj))
  dnf <- adj2dnf(adj*adj2)
  hierarchy <- list()
  vertices <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))
  children <- gsub(".*=", "", dnf)
  top <- hierarchy[[1]] <- vertices[-which(vertices %in% children)]
  vertices <- vertices[-which(vertices %in% top)]
  count <- 2
  while(length(vertices) > 0) {
    tmp <- dnf[grep(paste(paste(top, "=", sep = ""), collapse = "|"), dnf)]
    children <- unique(gsub(".*=", "", tmp))
    top <- hierarchy[[count]] <- children
    vertices <- vertices[-which(vertices %in% top)]
    count <- count + 1
  }
  return(hierarchy)
}

getHierarchyOld <- function(graph) {
  if (length(graph) == 0) {
    hierarchy <- NULL
  } else {
    adjmat <- dnf2adj(graph)
    get.order <- apply(adjmat, 2, sum)
    adjmat <- adjmat[order(get.order), order(get.order)]
    cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
    while(any(adjmat[lower.tri(adjmat)] == 1) & dim(cycles)[1] > 0) {
      for (i in 1:nrow(cycles)) {
        bad.cycles <- grep(paste(".*", rownames(adjmat)[cycles[i, 1]], ".*=", colnames(adjmat)[cycles[i, 2]], sep = ""), graph)
        if (length(bad.cycles) > 0) {
          graph <- graph[-bad.cycles]
          break()
        }
      }
      adjmat <- dnf2adj(graph)
      get.order <- apply(adjmat, 2, sum)
      adjmat <- adjmat[order(get.order), order(get.order)]
      cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
    }
    if (any(adjmat[lower.tri(adjmat)] == 1)) {
      cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
      for (i in 1:nrow(cycles)) {
        adjmat <- dnf2adj(graph)
        get.order <- apply(adjmat, 2, sum)
        adjmat <- adjmat[order(get.order), order(get.order)]
        cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1, arr.ind = TRUE)
        bad.cycles <- grep(paste(".*", rownames(adjmat)[cycles[i, 1]], ".*=", colnames(adjmat)[cycles[i, 2]], sep = ""), graph)
        if (length(bad.cycles) > 0) {
          graph <- graph[-bad.cycles]
        }
      }
    }
    nodes <- unique(unlist(strsplit(unlist(strsplit(graph, "=")), "\\+")))
    edges <- list()
    hierarchy <- list()
    children <- character()
    for (i in graph) {
      output <- unlist(strsplit(i, "="))
      edges[[output[2]]] <- c(edges[[output[2]]], unlist(strsplit(output[1], "\\+")))
      children <- unique(c(children, output[2]))
    }
    top <- nodes[-which(nodes %in% children)]
    left <- nodes
    lvl <- 1
    hierarchy[[lvl]] <- top
    while (length(left) > 0) {
      lvl <- lvl + 1
      left <- nodes[-which(nodes %in% unlist(hierarchy))]
      for (i in left) {
        if (all(edges[[i]] %in% unlist(hierarchy)) & !(i %in% unlist(hierarchy))) {
          if (length(hierarchy) < lvl) {
            hierarchy[[lvl]] <- i
          } else {
            if  (sum(edges[[i]] %in% hierarchy[[lvl]]) == 0) {
              hierarchy[[lvl]] <- c(hierarchy[[lvl]], i)
            }
          }
        }
      }
    }
  }
  for (i in 1:length(hierarchy)) {
    if (length(grep("!", hierarchy[[i]])) > 0) {
      hierarchy[[i]] <- hierarchy[[i]][-grep("!", hierarchy[[i]])]
    }
  }
  return(hierarchy)
}

getHierarchy3 <- function(graph) {
  if (length(graph) == 0) {
    hierarchy <- NULL
  } else {
    nodes <- unique(unlist(strsplit(unlist(strsplit(graph, "=")), "\\+")))
    edges <- list()
    hierarchy <- list()
    children <- character()
    for (i in graph) {
      output <- unlist(strsplit(i, "="))
      edges[[i]]$parents <- unlist(strsplit(output[1], "\\+"))
      edges[[i]]$child <- output[2]
      children <- unique(c(children, output[2]))
    }
    top <- nodes[-which(nodes %in% children)]
    left <- nodes
    lvl <- 1
    hierarchy[[lvl]] <- top
    while (length(left) > 0) {
      lvl <- lvl + 1
      left <- nodes[-which(nodes %in% unlist(hierarchy))]
      for (i in graph) {
        if (all(edges[[i]]$parents %in% unlist(hierarchy)) & !(edges[[i]]$child %in% unlist(hierarchy))) {
          if (length(hierarchy) < lvl) {
            hierarchy[[lvl]] <- edges[[i]]$child
          } else {
            if  (sum(edges[[i]]$parents %in% hierarchy[[lvl]]) == 0) {
              hierarchy[[lvl]] <- c(hierarchy[[lvl]], edges[[i]]$child)
            }
          }
        }
      }
    }
  }
  return(hierarchy)
}
  
getHierarchy2 <- function(graph) {
  adjmat <- dnf2adj(graph)
  if (is.null(adjmat)) {
    return(NULL)
  } else {
    get.order <- apply(abs(adjmat-1), 1, sum)
    get.order2 <- apply(adjmat, 2, sum)
    get.order3 <- get.order + get.order2
    adjmat <- adjmat[order(get.order2, decreasing = F), order(get.order, decreasing = F)]
    hierarchy <- list()
    lvl <- 1
    hierarchy[[lvl]] <- rownames(adjmat)[1]
    count <- 2
    nodes <- rownames(adjmat)
    nodes <- nodes[-which(nodes %in% rownames(adjmat)[1])]
    while(length(nodes) > 0) {
      if (all(adjmat[which(rownames(adjmat) %in% hierarchy[[lvl]]), ] == 0)) {
        hierarchy[[lvl]] <- c(hierarchy[[lvl]], rownames(adjmat)[count])
        nodes <- nodes[-which(nodes %in% rownames(adjmat)[count])]
      } else {
        lvl <- lvl + 1
        hierarchy[[lvl]] <- rownames(adjmat)[count]
        nodes <- nodes[-which(nodes %in% rownames(adjmat)[count])]
      }
      count <- count + 1
    }
    return(hierarchy)
  }
}
