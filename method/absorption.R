absorption <- function(bString, model) {
  graph <- model$reacID[which(bString == 1)]
  for (i in graph) {
    targets <- grep(paste("(?=.*", gsub("\\+", ")(?=.*", gsub("=", ")(?=.*=", i)), ")", sep = ""), graph, perl = TRUE)
    toomuch <- grep(paste("!", gsub("\\+", "|!", gsub("=.*", "", i)), "", sep = ""), graph[targets])
    if (length(toomuch) > 0) {
      targets <- targets[-grep(paste("!", gsub("\\+", "|!", gsub("=.*", "", i)), "", sep = ""), graph[targets])]
    }
    if (length(targets) > 1) {
      targets <- targets[-which(targets == which(graph %in% i))]
      bString[which(model$reacID %in% graph[targets])] <- 0
    }
  }
  return(bString)
}

absorptionII <- function(bString, model) {
  graph <- model$reacID[which(bString == 1)]
  nodes <- model$namesSpecies
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
      bString[which(model$reacID %in% graph[targets])] <- 0
    }
  }
  return(bString)
}

absorptionold <- function(bString, model) {
  graph <- model$reacID[bString == 1]
  gates <- list()
  for (i in graph) {
    if (is.null(gates[[i]])) {
      tmp <- unlist(strsplit(i, "="))
      output <- tmp[2]
      input <- unlist(strsplit(tmp[1], "\\+"))
      gates[[i]]$input <- input
      gates[[i]]$output <- output
    } else {
      input <- gates[[i]]$input
      output <- gates[[i]]$output
    }
    for (j in graph) {
      if (i %in% j) {
        next()
      }
      if (is.null(gates[[j]])) {
        tmp <- unlist(strsplit(j, "="))
        output2 <- tmp[2]
        input2 <- unlist(strsplit(tmp[1], "\\+"))
        gates[[j]]$input <- input2
        gates[[j]]$output <- output2
      } else {
        input2 <- gates[[j]]$input
        output2 <- gates[[j]]$output
      }
      if (all(input2 %in% input) & output %in% output2) {
        bString[which(model$reacID %in% i)] <- 0
      }
    }
  }
  return(bString)
}

absorptionIIold <- function(bString, model) {
  graph <- model$reacID[bString == 1]
  gates <- list()
  for (i in graph) {
    if (is.null(gates[[i]])) {
      tmp <- unlist(strsplit(i, "="))
      output <- tmp[2]
      input <- unlist(strsplit(tmp[1], "\\+"))
      gates[[i]]$input <- input
      gates[[i]]$output <- output
    } else {
      input <- gates[[i]]$input
      output <- gates[[i]]$output
    }
    for (j in graph) {
      if (i %in% j) {
        next()
      }
      if (is.null(gates[[j]])) {
        tmp <- unlist(strsplit(j, "="))
        output2 <- tmp[2]
        input2 <- unlist(strsplit(tmp[1], "\\+"))
        gates[[j]]$input <- input2
        gates[[j]]$output <- output2
      } else {
        input2 <- gates[[j]]$input
        output2 <- gates[[j]]$output
      }
      if (all(input %in% input2) & output %in% output2) {
        bString[which(model$reacID %in% i)] <- 0
      }
    }
  }
  return(bString)
}
