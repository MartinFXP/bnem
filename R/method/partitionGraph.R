partitionGraph <- function(model) { # does not yet work with cycles
  edges <- list()
  vertices <- list()
  count <- 0
  topnodes <- character()
  for (i in model$namesSpecies) {
    if (length(grep(paste("=", i, sep = ""), model$reacID)) == 0) {
      topnodes <- c(topnodes, i)
    }
  }
  botnodes <- character()
  for (i in model$namesSpecies) {
    if (length(grep(paste(".*", i, ".*=", sep = ""), model$reacID)) == 0) {
      botnodes <- c(botnodes, i)
    }
  }
  for (i in model$namesSpecies[-which(model$namesSpecies %in% topnodes)]) {
    count <- count + 1
    edges[[count]] <- model$reacID[grep(paste("=", i, sep = ""), model$reacID)]
    vertices[[count]] <- unique(unlist(strsplit(edges[[count]], "=")))
    vertices[[count]] <- gsub("!", "", unique(unlist(strsplit(vertices[[count]], "\\+"))))
  }
  ## try directed sequence and find which do not check out
  parents <- list()
  children <- list()
  for (i in 1:length(edges)) {
    parents[[i]] <- character()
    children[[i]] <- character()
    for (j in edges[[i]]) {
      tmp <- unlist(strsplit(j, "="))
      parents[[i]] <- c(parents[[i]], unlist(strsplit(tmp[1], "\\+")))
      children[[i]] <- c(children[[i]], tmp[2])
    }
    parents[[i]] <- unique(parents[[i]])
    children[[i]] <- unique(children[[i]])
  }
  ## detect cycles:
  ## graph <- model$reacID
  ## ddfs <- function(node, graph, path = NULL) {
  ##   graphCut <- graph[grep(paste("=", node, sep = ""), graph)]
  ##   for (i in graphCut) {
  ##     parents <- unlist(strsplit(i, "="))[1]
  ##     parents <- unlist(strsplit(parents, "\\+"))
  ##     for (j in parents) {
  ##       path <- c(path, node)
  ##       j <- gsub("!", "", j)
  ##       if (j %in% path) {
  ##         return(list(cycle = 1, path = path))
  ##       }
  ##       res <- ddfs(j, graph, path)
  ##       if (res$cycle == 1) {
  ##         return(list(cycle = 1, path = res$path))
  ##       } else {
  ##         cycle <- 0
  ##       }
  ##     }
  ##   }
  ##   return(list(cycle = 0, path = path))
  ## }
  ## cycles <- list()
  ## for (i in model$namesSpecies) {
  ##   tmp <- ddfs(i, graph)
  ##   if (tmp$cycle == 1) {
  ##     cycles[[i]] <- unique(tmp$path)
  ##   } else {
  ##     cycles[[i]] <- NULL
  ##   }
  ## }
  ## cycles.min <- cycles
  ## for (i in 1:length(cycles)) {
  ##   for (j in 1:length(cycles)) {
  ##     if (length(intersect(cycles[[i]],cycles[[j]])) == length(cycles[[i]]) & length(cycles[[i]]) != length(cycles[[j]])) {
  ##       cycles.min[[j]] <- NULL
  ##     }
  ##   }
  ## }
  ## order the stuff:
  
  suborder <- numeric()
  possible <- numeric(length(parents))
  while(min(possible) == 0) {
    for (i in 1:length(possible)) {
      if (min(possible[which(gsub("!", "", unlist(children)) %in% gsub("!", "", parents[[i]]))]) > 0) {
        possible[i] <- 1
        suborder <- c(suborder, i)
      }
    }
  }
  suborder <- unique(suborder)
  edges <- edges[suborder]
  vertices <- vertices[suborder]
  
  return(list(vertices = vertices, edges = edges))
}
