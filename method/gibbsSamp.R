gibbsSamp <- function(CNOlist=CNOlist, NEMlist=NEMlist, model=NCNOcutCompExp, approach=approach, parameters=parameters, verbose = TRUE, seeds = 1, iterations = 100, parallel=NULL) {
  if ((class(CNOlist) == "CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  if("abs" %in% approach) {
    if (length(NEMlist$norm) == 0) {
      print("data is not normed to (0,1); performing simple normalization")
      NEMlist$norm <- simpleNorm(NEMlist$exprs)
    }
  }
  if ("fc" %in% approach) {
    if (length(NEMlist$fc) == 0)  {
      print("foldchanges missing; automatic calculation")
      NEMlist$fc <- computeFc(CNOlist, NEMlist$exprs) 
      NEMlist <- computeSm(CNOlist, NEMlist, parameters)
    }
    if (length(NEMlist$E0) == 0)  {
      NEMlist <- computeSm(CNOlist, NEMlist, parameters)
    }
  }
  bLength <- length(model$reacID)
  simList = prep4sim(model)
  indexList = indexFinder(CNOlist, model)
  n = seeds # number starting strings

  if (!is.null(parallel)) {
    require(snowfall)
    #maxCores <- parallel::detectCores()
    #sfSetMaxCPUs(number=parallel)
    if (is.list(parallel)) {
      if (length(parallel[[1]]) != length(parallel[[2]])) { stop("The nodes (second list object in parallel) and the number of cores used on every node (first list object in parallel) must be the same length.") }
      n <- max(n, sum(parallel[[1]]))
      hosts <- character()
      for (i in 1:length(parallel[[1]])) {
        hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
      }
      hosts <- as.list(hosts)
      sfInit(parallel=TRUE, socketHosts=hosts)
    } else {
      n <- max(n, parallel)
      sfInit(parallel=TRUE, cpus=parallel)
    }
    sfExportAll()
    sfLibrary(CellNOptR)
    cat <- cat
    sfExport("cat")
  }
  # add a few random strings:
  bitStrings <- matrix(0, nrow = n, ncol = bLength+1)
  if (seeds >= 2) {
    for (i in 1:seeds) {
      if (round(i/2) == i/2) {
        bitStrings[i, ] <- abs(bitStrings[(i-1), ] - 1)
      } else {
        new <- FALSE
        safeNumber <- 0
        while(!new & safeNumber < 1000000000) {
          safeNumber <- safeNumber + 1
          bitStringTmp <- sample(c(0,1), bLength, replace = TRUE)
          compMatTmp <- t(t(bitStrings[, 1:bLength]) - bitStringTmp)
          compScore <- apply(abs(compMatTmp), 1, max)
          if (min(compScore) > 0) {
            bitStrings[i, 1:bLength] <- bitStringTmp
            new <- TRUE
          }
        }
      }
    }
  } else {
    bitStrings <- sample(c(0,1), bLength+1, replace = TRUE)
  }
  if (!is.matrix(bitStrings)) { bitStrings <- t(as.matrix(bitStrings)) }
  bitStringsSave <- numeric()
  checkSeed <- function(seed) {
    bitStringsSaveP <- numeric()
    startString <- bitStrings[seed, ]
    startString[length(startString)] <- computeScoreNemT1(CNOlist, model, startString[1:bLength], simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0)
    bitStringsSaveP <- rbind(bitStringsSaveP, startString)
    bitString <- startString
    for (i in 1:iterations) {
      timeMark <- Sys.time()
      for (j in 1:bLength) {
        typeZero <- bitString
        typeZero[j] <- 0
        typeOne <- bitString
        typeOne[j] <- 1
        typeZero[bLength+1] <- computeScoreNemT1(CNOlist, model, typeZero[1:bLength], simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0)
        typeOne[bLength+1] <- computeScoreNemT1(CNOlist, model, typeOne[1:bLength], simList, indexList, sizeFac, NAFac, approach = approach, NEMlist = NEMlist, parameters, tellme = 0)
        denomTmp <- typeZero[bLength+1] + typeOne[bLength+1]
        type <- sample(c(0,1), 1, prob = c((1 - typeZero[bLength+1]/denomTmp), (1 - typeOne[bLength+1]/denomTmp)))
        #bitStringsSaveP <- rbind(bitStringsSaveP, typeZero, typeOne)
        if (type == 0) {
          bitString <- typeZero
          bitStringsSaveP <- rbind(bitStringsSaveP, typeZero)
        } else {
          bitString <- typeOne
          bitStringsSaveP <- rbind(bitStringsSaveP, typeOne)
        }
      }
      timePassed <- as.numeric(Sys.time() - timeMark, unit = "secs")
      if (verbose) {
        print(paste("Iteration ", i, " (Seed: ", seed, ")", sep = ""))
        print(paste(" - Best String: ", toString(bitStringsSaveP[which(bitStringsSaveP[, ncol(bitStringsSaveP)] == min(bitStringsSaveP[, ncol(bitStringsSaveP)])[1]), ]), sep = ""))
        print(paste(" - Best Score: ", min(bitStringsSaveP[, ncol(bitStringsSaveP)]), sep = ""))
        print(paste(" - Iter_time: ", timePassed, sep = ""))
      }
    }
    return(bitStringsSaveP)
  }
  if (!is.null(parallel)) {
    seedResult <- sfLapply(as.list(1:seeds), checkSeed)
  } else {
    seedResult <- lapply(as.list(1:seeds), checkSeed)
  }
  for (i in 1:length(seedResult)) {
    bitStringsSave <- rbind(bitStringsSave, seedResult[[i]])
  }
  if (!is.null(parallel)) { sfStop() }
  # try to construct the best network from the sampled marginal distribution
  for (i in 1:bLength) {
  }
  return(bitStringsSave)
}
