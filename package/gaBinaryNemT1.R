gaBinaryNemT1 <- function (CNOlist,
                           model,
                           initBstring = NULL, # initBstring = TRUE
                           sizeFac = 1, 
                           NAFac = 1,
                           popSize = 50,
                           pMutation = 0.5,
                           maxTime = 60,
                           maxGens = 500, 
                           stallGenMax = 100,
                           relTol = 0.1, 
                           verbose = TRUE,
                           priorBitString = NULL,
                           selPress = c(1.2,0.0001), # 1.2
                           approach = "fc",
                           NEMlist,
                           fit = "linear",
                           targetBstring = "none",
                           elitism = NULL,
                           inversion = NULL,
                           graph = FALSE,
                           parameters = list(cutOffs = c(0.7,0.7,0.7), scoring = c(0.25,0.5,2)),
                           parallel = NULL, # parallel = N with N number of cores to use or a list with cores in the first and machines in the second entry like list(cores=c(2,4,8), machines=c("bionform1", "bioinform2", "bioinform3"))
                           parallel2 = 1,
                           selection = c("t"), # can be "t" or "s"
                           relFit = FALSE,
                           method = "none",
                           type = "SOCK",
                           exhaustive = FALSE,
                           delcyc = TRUE,
                           ...
                           ) {
  method <- checkMethod(method)
  if (parameters$cutOffs[1] > parameters$cutOffs[2]) {
    parameters$cutOffs <- sort(parameters$cutOffs)
    print(paste("your're cutoff parameters didn't make any sense. I can't let you do this, Dave. I changed them to ", parameters$cutOffs, ".", sep = ""))
  }
  if (is.null(elitism) == TRUE) { elitism <- 0 }
  if (elitism >= popSize) { elitism <- floor(0.1*popSize) }
  if (is.null(inversion) == TRUE) { inversion <- 0 }
  if (is.null(initBstring) == TRUE) {
    initBstring <- rep(1, length(model$reacID))
  }
  if (length(initBstring) == 1 & length(model$reacID) > 1) {
    initBstring <- max(0, min(initBstring, floor((popSize)*0.5)))
    initBstring <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist, model=model, approach=approach, seeds=initBstring, parameters=parameters, verbose = verbose, parallel=parallel, parallel2=parallel2,relFit = relFit,  method = method, ...)
    initBstring <- initBstring[order(localString[, ncol(localString)], decreasing = T), -ncol(initBstring)]
  }
  if ((class(CNOlist) == "CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  ## create foldchanges if not already included in nemlist
  CNOlist <- checkCNOlist(CNOlist)
  NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist, parameters = parameters, approach = approach, method)
  spaceExp <- 2^length(model$reacID)
  if ((popSize*stallGenMax) >= spaceExp & exhaustive) {
    print(paste("the genetic algorithm would score at least ", popSize*stallGenMax, " networks, while the size of the search space is only ", spaceExp, "; therefore an exhaustive search is initialised", sep = ""))
    result <- exSearch(CNOlist,model,sizeFac,NAFac,NEMlist,parameters,parallel=parallel,relFit = relFit,  method = method, ...)
    PopTolScores <- result$scores[which(result$scores < (result$score + abs(result$score)*relTol))]
    PopTol <- result$bStrings[which(result$scores < (result$score + abs(result$score)*relTol))]
    return(list(bString = result$bString, stringsTol = PopTol, stringsTolScores = PopTolScores, dtmRatio = result$dtmRatio, population = result$scores))
  } else {
    if ("s" %in% selection) {
      if (length(selPress) == 2) {
        selPressPct <- selPress[2]
        selPress <- selPress[1]
      } else {
        selPressPct <- 0
      }
      if (selPress < 1) {
        print("with selPress less than 1 low ranking networks are favoured")
        msg1 <- 1
      }
      if (selPress > 2 & fit == "linear") {
        selPress <- 2
        print("if selPress is greater than 2, fit cannot be set to linear; selPress set to 2; or restart with fit set to nonlinear")
        msg2 <- 1
      }
      if (selPress >= popSize) {
        selPress <- popSize - 1
        print(paste("selPress must be lower than popSize; selPress set to ", selPress, sep = ""))
        msg3 <- 1
      }
      if (selPress < 0) {
        selPress <- 1
        print("selPress has to be greater than zero and should be greater than 1; selPress set to 1")
        msg3 <- 1
      }
    }
    if ("t" %in% selection) {
      selPress <- floor(min(max(2, selPress[1]), popSize/2))
    }
    bLength <- length(model$reacID) # changed from length(initBstring)
    ## simList = prep4sim(model)
    simList <- NULL
    indexList = indexFinder(CNOlist, model)
    ## initialize starting population
    if (is.null(nrow(initBstring))) {
      initBstring <- t(as.matrix(initBstring))
    }
    Pop <- rbind(1 - initBstring, matrix(sample(c(0,1), (bLength * (popSize - (nrow(initBstring)*2))), replace = TRUE), nrow = (popSize - (nrow(initBstring)*2)), ncol = bLength), initBstring) # also the inverted initBstring is added
    Pop <- addPriorKnowledge(Pop, priorBitString)
    bestbit <- Pop[1, ]
    bestobj <- Inf
    stop <- FALSE
    g <- 0
    stallGen <- 0
    res <- rbind(c(g, bestobj, toString(bestbit), stallGen, Inf, 
                   Inf, toString(bestbit), 0), c(g, bestobj, toString(bestbit), 
                                                 stallGen, Inf, Inf, toString(bestbit), 0))
    colnames(res) <- c("Generation", "Best_score", "Best_bitString", 
                       "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen",
                       "Best_bit_Gen", "Iter_time")
    ## introduce the outputvector based on elitism on or off
    if (elitism >= 1) {
      res <- rbind(c(paste(res[4], " / ", res[1], sep = ""), paste(res[6], " (", res[5], ")", sep = ""),  res[7], res[8], "."))
      colnames(res) <- c("Stall Generations / Total Generations", "Best Score (Average Score)", "Best_bitString", "Iter_time", "----------------------------------")
    } else {
      res <- rbind(c(paste(res[4], " / ", res[1], sep = ""), paste(res[6], " (", res[5], ")", sep = ""),  res[7], res[2], res[3],res[8], "."))
      colnames(res) <- c("Stall Generations / Total Generations", "Best Score (Average Score)_Gen", "Best_bitString_Gen", "Best Score", "Best_bitString","Iter_time", "----------------------------------")
    }
    PopTol <- rep(NA, bLength)
    PopTolScores <- NA
    ## library(hash)
    ## scores2Hash = hash()
    ## getObj <- function(x) {
    ##   key = toString(x)
    ##   if (has.key(key, scores2Hash) == TRUE) {
    ##     return(scores2Hash[[key]])
    ##   } else {
    ##     Score = computeScoreNemT1(CNOlist, model, x, simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0)
    ##     if (length(scores2Hash) < 1000000) { # what's a reasonable number?
    ##       scores2Hash[[key]] = Score
    ##     }
    ##   }
    ##   return(Score)
    ## }
    getObj <- function(x) {
      Score <- computeScoreNemT1(CNOlist, model, x, simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, method = method, ...)
      return(Score)
    }
    t0 <- Sys.time()
    t <- t0
    selPressOrg <- selPress # needed to make the adjusted selection pressure work
    ## add a graph that shows the improvement:
    if (graph) {
      distances <- numeric()
      for (i in 1:100) {
        distances <- c(distances, dist(rbind(rep(1, i), rep(0, i))))
      }
      graphVal <- numeric()
      graphValAvg <- numeric()
      graphValSp <- numeric()
      gUp <- 0
      graphAxis <- numeric()
      graphics.off()
      par(bg = "white")
      graphValDist <- numeric()
      bestMemTurn <- 0
      bestMem <- "off"
    }
    if (!is.null(parallel)) {
      require(snowfall)
      if (is.list(parallel)) {
        if (length(parallel[[1]]) != length(parallel[[2]])) { stop("The nodes (second list object in parallel) and the number of cores used on every node (first list object in parallel) must be the same.") }
        hosts <- character()
        for (i in 1:length(parallel[[1]])) {
          hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
        }
        hosts <- as.list(hosts)
        sfInit(parallel=TRUE, socketHosts=hosts)
      } else {
        sfInit(parallel=TRUE, cpus=parallel, type = type)
      }
      sfLibrary(CellNOptR)
      sfExport(list = exportVars("ga"), local = T)
    }
    scores <- numeric(nrow(Pop))
    if ("t" %in% selection) {
      tsReduce <- function(x, scores) {
        y <- x[order(scores[x])[1]]
        return(y)
      }
    }
    if ("s" %in% selection) {
      susSel <- function(x, wheel1, breaks) {
        y <- which(wheel1 > breaks[x])[1]
        return(y)
      }
    }
    popMem <- numeric() # see high diversity later
    scoresMem <- numeric() # see high diversity later
    likelihoods <- numeric()
    lh.samples <- character()
    while (!stop) {
      if (delcyc) {
        if (!is.null(parallel)) {
          Pop <- t(sfApply(Pop, 1, removeCycles, model))
        } else {
          Pop <- t(apply(Pop, 1, removeCycles, model))
        }
      }
      bestReduced <- reduceGraph(Pop[popSize, ], model, CNOlist)
      if (sum((Pop[popSize, ] - bestReduced) != 0) > 0) {
        Pop[1, ] <- bestReduced
      }
      scores <- numeric(nrow(Pop))
      if (!is.null(parallel)) {
        scores <- sfApply(Pop, 1, getObj)
      } else {
        scores <- apply(Pop, 1, getObj)
      }
      scores[which(is.nan(scores) == TRUE)] <- max(scores[which(is.nan(scores) == FALSE)]) + 1 # a score could be NaN (probably because of negative feedback loops)
      rankP <- order(scores, decreasing = TRUE)
      Pop <- Pop[rankP, ]
      ## something to remember the samples:
      ## for (sample in 1:nrow(Pop)) {
      ##   if (!(toString(Pop[sample, ]) %in% lh.samples)) {
      ##     likelihoods <- c(likelihoods, scores[sample])
      ##     lh.samples <- c(lh.samples, toString(Pop[sample, ]))
      ##   }
      ## }
      ## try to alternatively keep in mind the best networks:
      if (graph) {
        if (bestMem == "off") {
          bestMem <- list()
          bestMem[[1]] <- Pop[popSize, ]
          bestMem[[2]] <- Pop[popSize, ]
          bestMemTurn <- abs(bestMemTurn - 1)
        } else {
          bestMem[[(bestMemTurn+1)]] <- Pop[popSize, ]
          bestMemTurn <- abs(bestMemTurn - 1)
        }
      }
      ## replace worst ranked networks with random networks for variation, ranking again would slow up the process and the idea is not to get random good networks but just good chunks in overall worse networks
      if (inversion >= 1) {
        Pop[1:inversion, ] <- abs(Pop[(popSize - inversion + 1):popSize, ] - 1) # replace the worst networks with the inverted best networks to create maximal diversity
      }
      scores <- scores[rankP]
      ## for high diversity, similar to sokolov & whitley (2005)
      ## popFull <- rbind(popMem, Pop)
      ## scoresFull <- c(scoresMem, scores)
      ## rankPFull <- order(scoresFull, decreasing = TRUE)
      ## popFull <- popFull[rankPFull, ]
      ## scoresFull <- scoresFull[rankPFull]
      ## bitsInDec <- apply(popFull, 1, bitToDec)
      ## uniquePop <- unique(bitsInDec)
      ## uniquePopPos <- numeric()
      ## for (i in length(scoresFull):(1 + max(popSize,(length(scoresFull) - length(uniquePop))))) {
      ##   uniquePopPos <- c(uniquePopPos, which(bitsInDec %in% uniquePop[i])[1])
      ## }
      ## uniquePopPos <- c(uniquePopPos, sample(1:nrow(popFull), (popSize - length(uniquePopPos))))
      ## popMem <- popFull[uniquePopPos, ]
      ## scoresMem <- scores[uniquePopPos]
      ## Pop <- popMem
      ## selection (s or t)
      PSize3 <- popSize - elitism
      if ("s" %in% selection) {
        ## nonlinear or linear fitness
        if (fit == "nonlinear") { # selPress can be higher than 2 in contrast to linear because in linear the chance of the worst member would get lower than 0 if the chance of the best would get higher than 2...
          x <- as.double(polyroot(c(rep(selPress, popSize-1), (selPress - popSize))))
          y <- polyroot(c(rep(selPress, popSize-1), (selPress - popSize)))
          x <- x[order(abs(y - as.double(y)))[1]]
          fitness <- (popSize*x^(1:popSize-1))/(sum(x^(1:popSize-1)))
        }
        if (fit == "linear") {
          fitness <- 2 - selPress + (2 * (selPress - 1) * (c(1:popSize) - 1)/(popSize - 1))
        }
        wheel1 <- cumsum(fitness/sum(fitness))
        breaks <- runif(1) * 1/popSize
        breaks <- c(breaks, breaks + ((1:(popSize - 1)))/popSize)
        sel <- rep(1, popSize)
        
        if (!is.null(parallel) & popSize > 10000) {
          sel <- sfApply(as.matrix(1:popSize), 1, susSel, wheel1, breaks)
        } else {
          sel <- apply(as.matrix(1:popSize), 1, susSel, wheel1, breaks)
        }
        ## for (i in 1:length(breaks)) {
        ##   sel[i] <- which(wheel1 > breaks[i])[1] # for a large population size, this could be parallelized
        ## }
      }
      if ("t" %in% selection) { # unbiased tournament selection sokolov & whitley (2005)
        
        pRanks <- sample(1:popSize, popSize)
        t.size <- min(popSize/2, selPress)
        ppRanks <- matrix(0, popSize, t.size)
        for (i in 1:t.size) {
          if (i == 1) {
            ppRanks[, i] <- pRanks[c(popSize, 1:(popSize-1))]
          } else {
            ppRanks[, i] <- ppRanks[c(popSize, 1:(popSize-1)), (i-1)]
          }
        }
        
        if (!is.null(parallel) & popSize > 10000) {
          sel <- as.vector(unlist(sfApply(cbind(pRanks, ppRanks), 1, tsReduce, scores)))
        } else {
          sel <- as.vector(unlist(apply(cbind(pRanks, ppRanks), 1, tsReduce, scores)))
        }
        
      }
      if ("r" %in% selection) {

        sel <- rep(1:popSize, round((1:popSize)/sum(1:popSize)*popSize))

        if (length(sel) < popSize) {
          sel <- c(sel, sel[length(sel)])
        }
        
      }
      if ("f" %in% selection) {

        scoresF <- scores*(-1)
        scoresF <- scoresF - min(scoresF)
        
        sel <- rep(1:popSize, round((scoresF)/sum(1:scoresF)*popSize))

        if (length(sel) < popSize) {
          sel <- c(sel, sel[1:(popSize-length(sel))])
        }
        
      }
      ##print(sel)
      ##print(scores)
      Pop2 <- Pop[sel, ]
      PSize2 <- dim(Pop2)[1]
      mates <- cbind(ceiling(runif(PSize3) * PSize2), ceiling(runif(PSize3) * 
                                                              PSize2))
      InhBit <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
                       ncol = bLength)
      InhBit <- InhBit < 0.5
      Pop3par1 <- Pop2[mates[, 1], ]
      Pop3par2 <- Pop2[mates[, 2], ]
      Pop3 <- Pop3par2
      Pop3[InhBit] <- Pop3par1[InhBit]
      MutProba <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
                         ncol = bLength)
      MutProba <- (MutProba < (pMutation/bLength))
      Pop3[MutProba] <- 1 - Pop3[MutProba]
      t <- c(t, Sys.time())
      g <- g + 1
      thisGenBest <- scores[popSize]
      thisGenBestBit <- Pop[popSize, ]
      if (is.na(thisGenBest)) {
        thisGenBest <- min(scores, na.rm = TRUE)
        thisGenBestBit <- Pop[which(scores == thisGenBest)[1], 
                              ]
      }
      if (thisGenBest < bestobj) {
        bestobj <- thisGenBest
        bestbit <- thisGenBestBit
        stallGen <- 0
      }
      else {
        stallGen <- stallGen + 1
      }
      ## try to get out of local minima by adjusting the selection pressure and introducing more diversity
      if ("s" %in% selection) {
        if (selPressPct < 0) {
          selPress <- max(selPress - selPress*selPressPct, 1)
        } else {
          if (fit == "linear") {
            selPress <- min(selPress + selPress*selPressPct, 2)
          } else {
            selPress <- min(selPress + selPress*selPressPct, popSize-1)
          }
        }
      }
      ## the following code needs changes if elitism is set to 0:
      resThisGen <- c(g, bestobj, toString(bestbit), stallGen, 
                      (mean(scores, na.rm = TRUE)), thisGenBest, toString(thisGenBestBit), 
                      as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"), "----------------------------------")
      names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", 
                             "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
                             "Best_bit_Gen", "Iter_time", "----------------------------------")
      ## verbose output based on elitism on or off
      if (elitism >= 1) {
        resThisGen <- c(g, bestobj, toString(bestbit), stallGen, 
                      (mean(scores, na.rm = TRUE)), thisGenBest, toString(thisGenBestBit), 
                      as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"), "----------------------------------")
        names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", 
                             "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
                             "Best_bit_Gen", "Iter_time", "----------------------------------")
        if (targetBstring[1] == "none") {
          resThisGen <- c(paste(resThisGen[4], " / ", resThisGen[1], sep = ""), paste(resThisGen[6], " (", resThisGen[5], ")", sep = ""),  resThisGen[7], resThisGen[8], "----------------------------------")
          names(resThisGen) <- c("Stall Generations / Total Generations", "Best Score (Average Score)", "Best_bitString", "Iter_time", "----------------------------------")
        } else {
          resThisGen <- c(paste(resThisGen[4], " / ", resThisGen[1], sep = ""), paste(resThisGen[6], " (", resThisGen[5], ")", sep = ""),  resThisGen[7], paste(targetBstring, collapse = ", "), resThisGen[8], "----------------------------------")
          names(resThisGen) <- c("Stall Generations / Total Generations", "Best Score (Average Score)", "Best_bitString", "Target_bitString", "Iter_time", "----------------------------------")
        }
      } else {
        resThisGen <- c(g, bestobj, toString(bestbit), stallGen, 
                      (mean(scores, na.rm = TRUE)), thisGenBest, toString(thisGenBestBit), 
                      as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"), "----------------------------------")
        names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", 
                             "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
                             "Best_bit_Gen", "Iter_time", "----------------------------------")
        if (targetBstring[1] == "none") {
          resThisGen <- c(paste(resThisGen[4], " / ", resThisGen[1], sep = ""), paste(resThisGen[6], " (", resThisGen[5], ")", sep = ""),  resThisGen[7], resThisGen[2], resThisGen[3],resThisGen[8], "----------------------------------")
          names(resThisGen) <- c("Stall Generations / Total Generations", "Best Score (Average Score)_Gen", "Best_bitString_Gen", "Best Score", "Best_bitString","Iter_time", "----------------------------------")
        } else {
          resThisGen <- c(paste(resThisGen[4], " / ", resThisGen[1], sep = ""), paste(resThisGen[6], " (", resThisGen[5], ")", sep = ""),  resThisGen[7], resThisGen[2], resThisGen[3], paste(targetBstring, collapse = ", "), resThisGen[8], "----------------------------------")
          names(resThisGen) <- c("Stall Generations / Total Generations", "Best Score (Average Score)_Gen", "Best_bitString_Gen", "Best Score", "Best_bitString", "Target_bitString","Iter_time", "----------------------------------")
        }
      }
      if (verbose) {
        if (elitism >= 1) {
          if (targetBstring[1] == "none") {
            print(paste(names(resThisGen)[3], ":", sep = ""))
            print(as.vector(resThisGen[3]))
            print(paste(names(resThisGen)[1], ": ", resThisGen[1], sep = ""))
            print(paste(names(resThisGen)[2], ": ", resThisGen[2], sep = ""))
            print(paste(names(resThisGen)[4], ": ", paste(resThisGen[4], "s @ ", Sys.time(), sep = ""), sep = ""))
            print(paste(names(resThisGen)[5], ": ", resThisGen[5], sep = ""))
          } else {
            print(paste(names(resThisGen)[3], ":", sep = ""))
            print(as.vector(resThisGen[3]))
            print(paste(names(resThisGen)[4], ":", sep = ""))
            print(as.vector(resThisGen[4]))
            print(paste(names(resThisGen)[1], ": ", resThisGen[1], sep = ""))
            print(paste(names(resThisGen)[2], ": ", resThisGen[2], sep = ""))
            print(paste(names(resThisGen)[5], ": ", resThisGen[5], sep = ""))
            print(paste(names(resThisGen)[6], ": ", paste(resThisGen[6], "s @ ", Sys.time(), sep = ""), sep = ""))
          }
        } else {
          print(resThisGen)
        }
      }
      res <- rbind(res, resThisGen)
      Criteria <- c((stallGen > stallGenMax), (as.numeric((t[length(t)] - 
                                                           t[1]), units = "secs") > maxTime), (g > maxGens))
      ## introduce a stop criteria for a target network (makes sense for simulation study)
      if (targetBstring[1] != "none" & g >= 2) {
        Criteria <- c((stallGen > stallGenMax), (as.numeric((t[length(t)] - 
                                                             t[1]), units = "secs") > maxTime), (g > maxGens), (computeScoreNemT1(CNOlist, model, bestbit, simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, ...) <= computeScoreNemT1(CNOlist, model, targetBstring, simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0, relFit = relFit, ...))) # put in here to stop if a network is found that has equal or better score than the target network
      }
      if (any(Criteria)) {
        stop <- TRUE
      }
      tolScore <- abs(scores[popSize] * relTol) # added the abs to allow for negative scores
      TolBs <- which(scores < scores[popSize] + tolScore)
      if (length(TolBs) > 0) {
        PopTol <- rbind(PopTol, Pop[TolBs, ])
        PopTolScores <- c(PopTolScores, scores[TolBs])
      }
      if (elitism > 0) {
        Pop <- rbind(Pop3, Pop[(popSize - elitism + 1):popSize, ])
      }
      else {
        Pop <- Pop3
      }
      Pop <- addPriorKnowledge(Pop, priorBitString)
      ## plot best network score, average networks score, selective pressure and network topology:
      if (graph) {
        split.screen(figs = c(1, 2), erase = FALSE) # screen 1 and 2 as columns
        split.screen(figs = c(2, 1), screen = 1, erase = FALSE) # screen 1 split in rows 3 and 4
        graphDraw <- 0
        if (length(graphVal) > 0) {
          if ((scores[order(scores, decreasing = TRUE)[popSize]] < graphVal[length(graphVal)]) | (stallGen %in% ceiling(stallGenMax*c(0.5,1)))) {
            #if (stallGen %in% ceiling(stallGenMax*c(0.5,1))) {
              ModelCut <- model
              ModelCut$interMat <- ModelCut$interMat[, as.logical(Pop[popSize, ])]
              ModelCut$notMat <- ModelCut$notMat[, as.logical(Pop[popSize, ])]
              ModelCut$reacID <- ModelCut$reacID[as.logical(Pop[popSize, ])]
              screen(2, new = T)
              plotDnf(ModelCut$reacID, CNOlist = CNOlist)
            #}
            graphVal <- c(graphVal, scores[order(scores, decreasing = TRUE)[popSize]])
            graphValAvg <- c(graphValAvg, sum(scores)/length(scores))
            graphValSp <- c(graphValSp, selPress)
            gUp <- gUp + 1
            graphAxis <- c(graphAxis, (g+1))
            graphDraw <- 1
            graphValDist <- c(graphValDist, dist(rbind(bestMem[[1]], bestMem[[2]])))
          }
        } else {
          graphVal <- c(graphVal, scores[order(scores, decreasing = TRUE)[popSize]])
          graphValAvg <- c(graphValAvg, sum(scores)/length(scores))
          graphValSp <- c(graphValSp, selPress)
          graphValDist <- c(graphValDist, dist(rbind(bestMem[[1]], bestMem[[2]])))
          gUp <- gUp + 1
          graphAxis <- c(graphAxis, (g+1))
          graphDraw <- 1
        }
        if (graphDraw == 1) {
          screen(3, new = T)
          plot(1:gUp, graphVal, col = "red", type = "l", main = paste("Score Improvement", sep = ""), xlab = "Generation", ylab = "Score", ylim = c(min(graphVal), max(c(graphVal, graphValAvg))), xaxt = "n")
          axis(1, at = 1:gUp, labels = graphAxis)
          legend(gUp, max(c(graphVal, graphValAvg)), legend = c("Best Score", "Average Score", "Selective Pressure"), fill = c("red", "blue", "black"), xjust = 1)
          lines(1:gUp, graphValAvg, col = "blue", type = "l")
          mtext("Selective Pressure", 4, line = 2)
          par(new=T)
          plot(1:gUp, graphValSp, col = "black", type = "l", axes=FALSE, ylab = "", xlab = "", yaxt = "n", ylim = c(1,max(2,ceiling(max(graphValSp)))))
          axis(4, at = (10:max(20,(ceiling(max(graphValSp))*10)))/10, labels = (10:max(20,(ceiling(max(graphValSp))*10)))/10)
          screen(4, new = T)
          plot(1:gUp, graphValDist, col = "black", type = "l", main = "Best Strings Distance", xlab = "Generation", ylab = "Distance", ylim = c(min(graphValDist), max(graphValDist)), xaxt='n')
          mtext("Edge Difference", 4, line = 2)
          abline(h = distances, lty = 3)
          axis(4, at = distances, labels = 1:length(distances))
          axis(1, at = 1:gUp, labels = graphAxis)
        }
      }
    }
    PopTol <- as.matrix(PopTol)[-1, ] # add the as.matrix to account for timeMax < "time needed for full generation"
    PopTolScores <- PopTolScores[-1]
    TolBs <- which(PopTolScores < scores[popSize] + tolScore)
    PopTol <- as.matrix(PopTol)[TolBs, ] # add the as.matrix to account for timeMax < "time needed for full generation"
    PopTolScores <- PopTolScores[TolBs]
    PopTolT <- cbind(PopTol, PopTolScores)
    PopTolT <- unique(PopTolT, MARGIN = 1)
    if (!is.null(dim(PopTolT))) {
      PopTol <- PopTolT[, 1:(dim(PopTolT)[2] - 1)]
      PopTolScores <- PopTolT[, dim(PopTolT)[2]]
    }
    else {
      PopTol <- PopTolT[1:(length(PopTolT) - 1)]
      PopTolScores <- PopTolT[length(PopTolT)]
    }
    res <- res[2:dim(res)[1], ]
    rownames(res) <- NULL
    if (!is.null(parallel)) {
      sfStop()
    }
    bestbit <- reduceGraph(bestbit, model, CNOlist)
    dtmRatio <- 0#mean(abs(likelihoods - mean(likelihoods)))/abs(min(likelihoods) - mean(likelihoods))
    #names(likelihoods) <- lh.samples
    return(list(bString = bestbit, results = res, stringsTol = PopTol, stringsTolScores = PopTolScores, dtmRatio = dtmRatio, population = likelihoods))
  }
}
