#' simulates the activation pattern (truth table) of a hyper-graph and annotated perturbation experiments
#' @param CNOlist, CNOlist object
#' @param model model object
#' @param bString binary vector denoting the sub-graph given model
#' @param NEMlist NEMlist object only for devel
#' @author Martin Pirkl
#' @return return the truth tables for certain perturbation experiments
#' @export
#' @import
#' CellNOptR
#' matrixStats
#' @examples
#' library(bnem)
#' library(CellNOptR)
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
#' model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' states <- simulateStatesRecursive(CNOlist, model, rep(1, length(model$reacID)))
simulateStatesRecursive <-
function(CNOlist, model, bString, NEMlist = NULL) {
  getState <- function(CNOlist, node, signalStates, graph, children = NULL, NEMlist = NULL) {
    graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
    if (length(graphCut) == 0) {
      signalStates[, node] <- 0
    } else {
      sop <- numeric(nrow(signalStates))
      children2 <- gsub("!", "", children)
      for (i in graphCut) {
        parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
        if (length(parents) == 0) {
          pob <- rep(0, nrow(signalStates))
        } else {
          pob <- rep(1, nrow(signalStates))
        }
        for (j in parents) {
          j2 <- gsub("!", "", j)
          if (any(is.na(signalStates[, j2]) == T)) {
            if (j %in% j2) {
              node2 <- node
              add1 <- 0
            } else {
              node2 <- paste("!", node, sep = "")
              add1 <- 1
            }
            if (j2 %in% children2) {
              ## this speeds up the process and will
              subGraph <- graph
              subGraph <- subGraph[-grep(paste(".*=", node, "|.*", j2, ".*=.*", sep = ""), subGraph)]
              signalStatesTmp <- getState(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = subGraph, children = children2[-which(children2 %in% node)], NEMlist)
              
              ## ## is mathematically more reasonable but takes longer:
              ## signalStates2 <- signalStates
              ## signalStates2[, node] <- 0
              ## signalStatesTmp <- getState(CNOlist = CNOlist, node = j2, signalStates = signalStates2, graph = subGraph, children = children2[-which(children2 %in% node)], NEMlist)
              
              ## if ((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2 != ceiling((length(grep("!", children[which(children2 %in% j2):length(children2)]))+add1)/2)) {
              ##   ## negative feedback loop calculation does not seem to be general enough and also not feasible:
              ## } else {
              ## }
              if (add1 == 0) {
                pobMult <- signalStatesTmp[, j2]
              } else {
                pobMult <- add1 - signalStatesTmp[, j2]
              }
            } else {
              signalStates <- getState(CNOlist = CNOlist, node = j2, signalStates = signalStates, graph = graph, children = unique(c(children, node2)), NEMlist)
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
      if (node %in% colnames(CNOlist@inhibitors)) {
        sop <- sop*(1 - CNOlist@inhibitors[, node])
      }
      if (node %in% colnames(CNOlist@stimuli)) {
        sop <- max(sop, CNOlist@stimuli[, node])
      }
      signalStates[, node] <- sop
    }
    return(signalStates)
  }
  bString <- reduceGraph(bString, model, CNOlist)
  stimuli <- colnames(CNOlist@stimuli)
  signals <- sort(c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))]))
  graph0 <- model$reacID[which(bString == 1)]
  stimuliStates <- CNOlist@stimuli
  if (!is.null(NEMlist$signalStates)) {
    signalStates <- NEMlist$signalStates
  } else {
    signalStates <- matrix(NA, nrow = nrow(CNOlist@signals[[2]]), ncol = length(signals))
    rownames(signalStates) <- rownames(CNOlist@signals[[2]])
    colnames(signalStates) <- signals
    signalStates <- cbind(stimuliStates, signalStates)
  }
  for (k in signals) {
    if (sum(is.na(signalStates[, k]) == T) == length(signalStates[, k])) {
      signalStates <- getState(CNOlist = CNOlist, node = k, signalStates = signalStates, graph = graph0, children = NULL, NEMlist)
    }
  }
  signalStates <- signalStates[, which(colnames(signalStates) %in% colnames(CNOlist@signals[[1]]))]
  if (ncol(CNOlist@signals[[1]]) != 1) {
    signalStates <- signalStates[, order(colnames(signalStates))]
  } else {
    signalStates <- as.matrix(signalStates)
    colnames(signalStates) <- colnames(CNOlist@signals[[1]])
  }
  return(signalStates = signalStates)
}
