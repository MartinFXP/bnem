#' @noRd
#' @import matrixStats
simulateDnf <-
function(dnf, stimuli = NULL, inhibitors = NULL) {
  getStateDnf <- function(node, signalStates, graph, children = NULL) {
    graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
    if (length(graphCut) == 0) {
      signalStates[, node] <- 0
    } else {
      sop <- numeric(nrow(signalStates))
      children2 <- gsub("!", "", children)
      for (i in graphCut) {
        parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
        pob <- rep(1, nrow(signalStates))
        for (j in parents) {
          j2 <- gsub("!", "", j)
          if (sum(is.na(signalStates[, j2]) == T) == length(signalStates[, j2])) {
            if (j %in% j2) {
              node2 <- node
              add1 <- 0
            } else {
              node2 <- paste("!", node, sep = "")
              add1 <- 1
            }
            if (j2 %in% children2) {
                                        #subGraph <- graph[-grep(paste(".*", j2, ".*=", node, sep = ""), graph)]
                                        #subGraph <- graph[-grep(paste(".*=", node, sep = ""), graph)] # seems to work
              subGraph <- graph[-grep(paste(".*=", node, "|.*", j2, ".*=.*", sep = ""), graph)] # not needed is it? I think it is!!!
              signalStatesTmp <- getStateDnf(node = j2, signalStates = signalStates, graph = subGraph, children = NULL)
              if ((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2 != ceiling((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2)) {
                ## negative feedback loop calculation does not seem to be general enough and also not feasible:
              } else {
              }
              if (add1 == 0) {
                pobMult <- signalStatesTmp[, j2]
              } else {
                pobMult <- add1 - signalStatesTmp[, j2]
              }
            } else {
              signalStates <- getStateDnf(node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)))
              if (add1 == 0) {
                pobMult <- signalStates[, j2]
              } else {
                pobMult <- add1 - signalStates[, j2]
              }
            }
            pob <- pob*pobMult
          } else {
            if (j %in% j2) {
              add1 <- 0
            } else {
              add1 <- 1
            }
            if (add1 == 0) {
              pobMult <- signalStates[, j2]
            } else {
              pobMult <- add1 - signalStates[, j2]
            }
            pob <- pob*pobMult
          }
          if (max(pob, na.rm = T) == 0) { break() }
        }
        sop <- sop + pob
        if (min(sop, na.rm = T) > 0) { break() }
      }
      sop[sop > 0] <- 1
      if (node %in% inhibitors) {
        sop <- sop*0
      }
      if (node %in% stimuli) {
        sop <- max(sop, 1)
      }
      signalStates[, node] <- sop
    }
    return(signalStates)
  }
  signals <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))
  graph <- dnf
  signalStates <- matrix(NA, nrow = 1, ncol = length(signals))
  rownames(signalStates) <- paste(c("stimuli:", stimuli, "inhibitors:", inhibitors), collapse = " ")
  colnames(signalStates) <- signals
  signalStates[which(signals %in% stimuli)] <- 1
  for (k in signals) {
    if (is.na(signalStates[, k]) == T) {
      signalStates <- getStateDnf(node = k, signalStates = signalStates, graph = graph, children = NULL)
    }
  }
  namestmp <- colnames(signalStates)
  signalStates <- as.vector(signalStates)
  names(signalStates) <- namestmp
  return(signalStates = signalStates)
}
