localSearch <- function(CNOlist, NEMlist, model, approach = "fc", initSeed = NULL, seeds = 1, parameters = list(cutOffs = c(0,1,0), scoring = c(0.25,0.5,2)), sizeFac = 10^-10, NAFac = 1, relTol = 0.01, verbose = TRUE, parallel=NULL, parallel2 = 1, relFit = FALSE, method = "s", max.steps = Inf, max.time = Inf, node = NULL, absorpII = TRUE, draw = TRUE, prior = NULL) {
  require(matrixStats)
  if (is.null(prior)) {
    prior <- rep(0, length(model$reacID))
  }
  method <- checkMethod(method)
  debug <- F
  if (verbose %in% "part") {
    verbose2 <- "part"
    verbose <- TRUE
  }
  if (seeds == "max") {
    seeds <- length(model$reacID)*2
    seeds2 <- "max"
  } else {
    seeds2 <- "none"
  }
  if ((class(CNOlist) == "CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  CNOlist <- checkCNOlist(CNOlist)
  NEMlist <- checkNEMlist(NEMlist, CNOlist, parameters, approach, method)
  bLength <- length(model$reacID)
  ##simList = prep4sim(model)
  indexList = indexFinder(CNOlist, model)
  n = seeds # number of strings to analyse
  bitStrings <- matrix(0, nrow = n, ncol = bLength)
  if (n >= 1) {
    bitStrings[1, ] <- c(0, rep(0, bLength-1))
  }
  if (!is.null(initSeed)) {
    bitStrings[1, ] <- reduceGraph(initSeed, model, CNOlist)
  }
  ## add a few random (but good) strings:
  if (seeds >= 2 & seeds2 != "max") {
    bitStrings[2:seeds, ] <- createCube(ncol(bitStrings), seeds-1) # creates network far opposite each other!
  }
  ## if (seeds >= 2 & seeds2 != "max") {
  ##   for (i in 2:seeds) {
  ##     if (round(i/2) == i/2) {
  ##       bitStrings[i, ] <- reduceGraph(abs(bitStrings[(i-1), ] - 1), model, CNOlist)
  ##     } else {
  ##       new <- FALSE
  ##       safeNumber <- 0
  ##       while(!new & safeNumber < 1000000000) {
  ##         safeNumber <- safeNumber + 1
  ##         bitStringTmp <- reduceGraph(sample(c(0,1), bLength, replace = TRUE), model, CNOlist)
  ##         compMatTmp <- t(t(bitStrings) - bitStringTmp)
  ##         ##compScore <- apply(abs(compMatTmp), 1, max)
  ##         compScore <- rowMaxs(abs(compMatTmp))
  ##         if (min(compScore) > 0) {
  ##           bitStrings[i, ] <- bitStringTmp
  ##           new <- TRUE
  ##         }
  ##       }
  ##     }
  ##   }
  ## }
  ## if (seeds >= 2 & seeds2 == "max") {
  ##   for (i in 2:seeds) {
  ##     if (round(i/2) == i/2) {
  ##       bitStrings[i, ] <- abs(bitStrings[(i-1), ] - 1)
  ##     } else {
  ##       if (i == 1) {
  ##         bitStrings[i, ] <- sample(c(0,1), length(model$reacID), replace = T)
  ##       }
  ##       for (j in 1:(i-1)) {
  ##         bitStrings[i, ] <- bitStrings[i-2, ((j-1)*(length(model$reacID)/(i-1)) + 1):(j*(length(model$reacID)/(i-1)))]
  ##         #bitStrings[i, ] <- c(bitStrings[i-2, 1:(seeds/(i-1))], bitStrings[i-1, 1:(seeds/(i-1))])
  ##       }
  ##     }
  ##   }
  ## }
  bitStringsMem <- numeric()
  bitStringsScores <- numeric()
  if (!is.null(parallel)) {
    require(snowfall)
                                        #maxCores <- parallel::detectCores()
                                        #sfSetMaxCPUs(number=parallel)
    if (is.list(parallel)) {
      if (length(parallel[[1]]) != length(parallel[[2]])) { stop("The nodes (second list object in parallel) and the number of cores used on every node (first list object in parallel) must be the same.") }
      hosts <- character()
      for (i in 1:length(parallel[[1]])) {
        hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
      }
      hosts <- as.list(hosts)
      sfInit(parallel=TRUE, cpus=sum(parallel[[1]]), type="SOCK", socketHosts=hosts)
    } else {
      sfInit(parallel=TRUE, cpus=parallel)
    }
    sfLibrary(CellNOptR)
    sfExport(list = exportVars("loc"), local = T)
  }
  start.time <- Sys.time()
  processSeed <- function(i, bitStrings, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method, max.steps = max.steps, node = node) {
    stimuli <- colnames(CNOlist@stimuli)
    inhibitors <- colnames(CNOlist@inhibitors)
    signals <- colnames(CNOlist@signals[[1]])
    step.count <- 0
    row <- i
    new <- FALSE
    bitString <- bitStrings[row, ]
    compMatTmp <- t(t(bitStringsMem) - bitString)
    ##compScore <- apply(abs(compMatTmp), 1, max)
    compScore <-rowMaxs(abs(compMatTmp))
    if (min(compScore) == 0) {
      new <- FALSE
      safeNumber <- 0
      while(!new & safeNumber < 1000) {
        safeNumber <- safeNumber + 1
        bitString <- sample(c(0,1), bLength, replace = TRUE)
        compMatTmp <- t(t(rbind(bitStrings,bitStringsMem)) - bitStringTmp)
        ##compScore <- apply(abs(compMatTmp), 1, max)
        compScore <- rowMaxs(abs(compMatTmp))
        if (min(compScore) > 0) {
          new <- TRUE
        }
      }
      print(safeNumber)
    }
    fullScore <- computeScoreNemT1(CNOlist, model = model, bitString, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
    if (verbose) {
      print(paste("Seed Network ", row, "/", n, sep = ""))
      if (!(verbose %in% "part")) {
        print(toString(bitString))
      }
      print(paste(" - Score: ", fullScore, sep = ""))
      print("--------------------------------------------------")
      if (any(bitString != 0) & draw) {
        plotDnf(model$reacID[which(bitString == 1)], CNOlist = CNOlist)
        ##plotBinary(bitString, model, signals = signals, inhibitors = inhibitors, stimuli = stimuli)
      }
    }
    ## put edges that were just changed (added or deleted) on a "tabu"-list to prevent them from being changed in the next n iterations. atm: no.
    stop <- FALSE #
    save.scores <- numeric()
    edges.changed <- numeric()
    edge.history <- character()
    counter <- 0
    while(!stop) {
      score <- computeScoreNemT1(CNOlist, model = model, bitString, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
      save.scores <-c(save.scores, score)
      scores <- numeric(bLength)
      sizes <- numeric(bLength)
      pValues <- numeric(bLength) + 1
      timeMark <- Sys.time()
      scoreThem <- function(i, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method) {
        if (!is.null(node)) {
          if (length(grep(paste(node, collapse = "|"), model$reacID[i])) == 0) {
            return(Inf)
          } else {
            bitStringTmp <- bitString
            bitStringTmp[i] <- abs(1 - bitStringTmp[i])
            redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
            if (all(bitString == redBstring) & absorpII) {
              bitStringTmp <- absorptionII(bitStringTmp, model)
            } else {
              bitStringTmp <- redBstring
            }
            return(computeScoreNemT1(CNOlist, model = model, bitStringTmp, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method))
          }
        } else {
          bitStringTmp <- bitString
          bitStringTmp[i] <- abs(1 - bitStringTmp[i])
          redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
          if (all(bitString == redBstring) & absorpII) {
            bitStringTmp <- absorptionII(bitStringTmp, model)
          } else {
            bitStringTmp <- redBstring
          }
          if (sizeFac == 0) {
            size <- length(unlist(strsplit(model$reacID[which(bitStringTmp == 1)], "\\+")))
            return(c(computeScoreNemT1(CNOlist, model = model, bitStringTmp, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method), size))
          } else {
            return(computeScoreNemT1(CNOlist, model = model, bitStringTmp, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method))
          }
        }
      }
      edge.matrix <- as.matrix(1:bLength)
      if (sum(prior != 0) > 0) {
        edge.matrix <- as.matrix(edge.matrix[-which(prior != 0), ])
      }
      if (parallel2 == 1 & !is.null(parallel)) {
        scores <- sfApply(edge.matrix, 1, scoreThem, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method)
      } else {
        scores <- apply(edge.matrix, 1, scoreThem, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method)
      }
      if (sizeFac == 0) {
        sizes <- scores[2, ]
        scores <- scores[1, ] 
        size <- length(unlist(strsplit(model$reacID[which(bitString == 1)], "\\+")))
        check.size <- FALSE
        if (any(scores == score) & all(scores >= score)) {
          if (any(sizes[which(scores == score)] < size)) {
            check.size <- TRUE
          }
        }
      } else {
        check.size <- FALSE
      }
      if (sum(prior != 0) > 0) {
        scores.tmp <- scores
        scores <- numeric(length(model$reacID))
        scores[which(prior == 0)] <- scores.tmp
        scores[which(prior != 0)] <- Inf
      }
      ##   for (i in 1:bLength) {
      ##     if (debug) {
      ##       print(model$reacID[i])
      ##     }
      ##     if (!is.null(node)) {
      ##       if (length(grep(paste(node, collapse = "|"), model$reacID[i])) == 0) {
      ##         scores[i] <- Inf
      ##       } else {
      ##         bitStringTmp <- bitString
      ##         bitStringTmp[i] <- abs(1 - bitStringTmp[i])
      ##         redBstring <- reduceGraph(bitStringTmp, model)
      ##         if (all(bitString == redBstring) & absorpII) {
      ##           bitStringTmp <- absorptionII(bitStringTmp, model)
      ##         } else {
      ##           bitStringTmp <- redBstring
      ##         }
      ##         scores[i] <- computeScoreNemT1(CNOlist, model = model, bitStringTmp, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
      ##       }
      ##     } else {
      ##       bitStringTmp <- bitString
      ##       bitStringTmp[i] <- abs(1 - bitStringTmp[i])
      ##       redBstring <- reduceGraph(bitStringTmp, model)
      ##       if (all(bitString == redBstring) & absorpII) {
      ##         bitStringTmp <- absorptionII(bitStringTmp, model)
      ##       } else {
      ##         bitStringTmp <- redBstring
      ##       }
      ##       scores[i] <- computeScoreNemT1(CNOlist, model = model, bitStringTmp, sizeFac = sizeFac, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, method = method)
      ##     }
      ##   }
      ## }
      
      ## data <- scores - score
      ## for (i in colnames(CNOlist@inhibitors)) {
      ##   toGrep <- paste("=", i, sep = "")
      ##   toGrep <- i
      ##   dataGrep <- data[grep(toGrep, model$reacID)]
      ##   for (j in grep(toGrep, model$reacID)) {
      ##     tValueTmp <- (mean(dataGrep) - data[j]) / (sd(dataGrep)/sqrt(length(dataGrep)))
      ##     pValues[j] <- (1 - pt(tValueTmp, df=length(dataGrep) - 1))
      ##   }
      ## }
      timePassed <- as.numeric(Sys.time() - timeMark, unit = "secs")
      if (sum(scores < score) > 0 | check.size) {
        if (check.size) { ##pValues[is.na(pValues)] <- 1
          pValues <- scores # if you do not want to use pValues
          topGates <- which(scores == score)
          topScores <- scores[topGates]
          topPvalues <- pValues[topGates]
          minPvalue <- min(topPvalues)
          topGate <- which(topPvalues == minPvalue)
          topGate <- topGate[which(sizes[topGate] == min(sizes[topGate]))[1]]
          topScore <- topScores[topGate]
          whichGate <- which(sizes < size & scores == score)[1] # topGates[topGate]
          bitStringTmp <- bitString
          bitStringTmp[whichGate] <- abs(1 - bitStringTmp[whichGate])
          redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
          if (all(bitString == redBstring) & absorpII) {
            bitStringTmp2 <- bitString
            bitString <- absorptionII(bitStringTmp, model)
            edges <- sum(bitString != bitStringTmp2)
            edges.changed <- c(edges.changed, -edges)
            if (verbose) {
              deleted <- which(bitString == 0 & bitStringTmp == 1)
              if (bitStringTmp[whichGate] == 0 & length(deleted) > 1) {
                deleted <- deleted[-which(deleted == whichGate)]
              }
              print(paste("Edges changed: ", edges, sep = ""))
              print(paste("Deleted gates due to inverse absorption: ", paste(model$reacID[deleted], collapse = ", "), sep = ""))
            }
          } else {
            bitStringTmp2 <- bitString
            bitString <- redBstring
            deleted <- which(bitString == 0 & bitStringTmp2 == 1)
            edges <- sum(bitString != bitStringTmp2)
            edges.changed <- c(edges.changed, edges)
            if (verbose & (length(deleted) > 1 | (length(deleted) > 0 & bitString[whichGate] == 1))) {
              print(paste("Edges changed: ", edges, sep = ""))
              print(paste("Deleted gates due to absorption: ", paste(model$reacID[deleted], collapse = ", "), sep = ""))
            }
          }
          bitStringsMem <- rbind(bitStringsMem, bitString)
          if (bitString[whichGate] == 1) {
            edge.history <- c(edge.history, model$reacID[whichGate])
          }
          if (verbose) {
            counter <- counter + 1
            print(paste("Iter step: ", counter, sep = ""))
            if (bitString[whichGate] == 1) {
              print(paste("Added gate ", model$reacID[whichGate], sep = ""))
              if (!(verbose %in% "part")) {
                print(toString(bitString))
              }
              print(paste(" - Score: ", topScore, sep = ""))
              print(paste(" - Iter_time: ", timePassed, " @ ", Sys.time(), sep = ""))
              print("--------------------------------------------------")
            } else {
              print(paste("Deleted gate ", model$reacID[whichGate], sep = ""))
              if (!(verbose %in% "part")) {
                print(toString(bitString))
              }
              print(paste(" - Score: ", topScore, sep = ""))
              print(paste(" - Iter_time: ", timePassed, " @ ", Sys.time(), sep = ""))
              print("--------------------------------------------------")
            }
            if (any(bitString != 0) & draw) {
              plotDnf(model$reacID[which(bitString == 1)], CNOlist = CNOlist)
            }
          }
          step.count <- step.count + 1
          if (step.count >= max.steps) {
            if (verbose) print("no more steps")
            stop <- TRUE
          }
        } else {
          ##pValues[is.na(pValues)] <- 1
          pValues <- scores # if you do not want to use pValues
          topGates <- which(scores < score)
          topScores <- scores[topGates]
          topPvalues <- pValues[topGates]
          minPvalue <- min(topPvalues)
          topGate <- which(topPvalues == minPvalue)[1]
          topScore <- topScores[topGate]
          whichGate <- topGates[topGate]
          bitStringTmp <- bitString
          bitStringTmp[whichGate] <- abs(1 - bitStringTmp[whichGate])
          redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
          if (all(bitString == redBstring) & absorpII) {
            bitStringTmp2 <- bitString
            bitString <- absorptionII(bitStringTmp, model)
            edges <- sum(bitString != bitStringTmp2)
            edges.changed <- c(edges.changed, -edges)
            if (verbose) {
              deleted <- which(bitString == 0 & bitStringTmp == 1)
              if (bitStringTmp[whichGate] == 0 & length(deleted) > 1) {
                deleted <- deleted[-which(deleted == whichGate)]
              }
              print(paste("Edges changed: ", edges, sep = ""))
              print(paste("Deleted gates due to inverse absorption: ", paste(model$reacID[deleted], collapse = ", "), sep = ""))
            }
          } else {
            bitStringTmp2 <- bitString
            bitString <- redBstring
            deleted <- which(bitString == 0 & bitStringTmp2 == 1)
            edges <- sum(bitString != bitStringTmp2)
            edges.changed <- c(edges.changed, edges)
            if (verbose & (length(deleted) > 1 | (length(deleted) > 0 & bitString[whichGate] == 1))) {
              print(paste("Edges changed: ", edges, sep = ""))
              print(paste("Deleted gates due to absorption: ", paste(model$reacID[deleted], collapse = ", "), sep = ""))#paste(model$reacID[deleted[-which(deleted == whichGate)]], collapse = ", "), sep = ""))
            }
          }
          bitStringsMem <- rbind(bitStringsMem, bitString)
          if (bitString[whichGate] == 1) {
            edge.history <- c(edge.history, model$reacID[whichGate])
          }
          if (verbose) {
            counter <- counter + 1
            print(paste("Iter step: ", counter, sep = ""))
            if (bitString[whichGate] == 1) {
              print(paste("Added gate ", model$reacID[whichGate], sep = ""))
              if (!(verbose %in% "part")) {
                print(toString(bitString))
              }
              print(paste(" - Score: ", topScore, sep = ""))
                                        #print(paste(" - P.Value: ", minPvalue, sep = ""))
              print(paste(" - Iter_time: ", timePassed, " @ ", Sys.time(), sep = ""))
              print("--------------------------------------------------")
            } else {
              print(paste("Deleted gate ", model$reacID[whichGate], sep = ""))
              if (!(verbose %in% "part")) {
                print(toString(bitString))
              }
              print(paste(" - Score: ", topScore, sep = ""))
                                        #print(paste(" - P.Value: ", minPvalue, sep = ""))
              print(paste(" - Iter_time: ", timePassed, " @ ", Sys.time(), sep = ""))
              print("--------------------------------------------------")
            }
            if (any(bitString != 0) & draw) {
              plotDnf(model$reacID[which(bitString == 1)], CNOlist = CNOlist)
              ##plotBinary(bitString, model, signals = signals, inhibitors = inhibitors, stimuli = stimuli)
            }
          }
          step.count <- step.count + 1
          if (step.count >= max.steps) {
            if (verbose) print("no more steps")
            stop <- TRUE
          }
        }
      } else {
        if (verbose) print("no further improvement")
        stop <- TRUE
      }
      if (Sys.time() - start.time > max.time) {
        if (verbose) print("no more time")
        stop <- TRUE
      }
    }
    return(list(A=score,B=bitString,row=row, C=save.scores, D=edges.changed, E=edge.history))
  }
  if (!is.null(parallel) & parallel2 == 0) {
    bsTemp <- sfApply(as.matrix(1:n), 1, processSeed, bitStrings, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method, max.steps = max.steps, node = node)
  } else {
    bsTemp <- apply(as.matrix(1:n), 1, processSeed, bitStrings, CNOlist, model, simList = NULL, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method, max.steps = max.steps, node = node)
  }
  save.scores <- list()
  edges.changed <- list()
  edge.history <- list()
  for (i in 1:length(bsTemp)) {
    bitStringsScores[bsTemp[[i]]$row] <- bsTemp[[i]]$A
    bitStrings[bsTemp[[i]]$row, ] <- bsTemp[[i]]$B
    save.scores[[i]] <- bsTemp[[i]]$C
    edges.changed[[i]] <- bsTemp[[i]]$D
    edge.history[[i]] <- bsTemp[[i]]$E
  }
  ## bitString <- bitStrings[which(bitStringsScores == min(bitStringsScores))[1], ] # for single best bitString
  bitString <- bitStrings#[which(bitStringsScores <= (min(bitStringsScores) + abs(min(bitStringsScores)*relTol))), ]
  if (!is.null(parallel)) {
    sfStop()
  }
  bitString <- bitStrings # cbind(bitStrings, bitStringsScores)
  return(list(bStrings = bitString, scores = save.scores, edges = edges.changed, history = edge.history))
}
