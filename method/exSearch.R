exSearch <- function(CNOlist,model,sizeFac=10^-10,NAFac=1,NEMlist,parameters=list(cutOffs=c(0,1,0), scoring=c(0.1,0.2,0.9)), parallel = NULL, method = "s", relFit = F, verbose = TRUE, reduce = TRUE, ...) {

  CNOlist <- checkCNOlist(CNOlist)
  
  cutModel2 <- function (model, bString) {
    if (sum(bString == 1) > 0) {
      bs = as.logical(bString)
      newmodel <- list()
      if (is.null(dim(model$interMat))) {
        newmodel$interMat <- model$interMat[bs]
        newmodel$notMat <- model$notMat[bs]
        newmodel$reacID <- model$reacID[bs]
        newmodel$namesSpecies <- model$namesSpecies
      } else {
        newmodel$interMat <- model$interMat[, bs]
        newmodel$notMat <- model$notMat[, bs]
        newmodel$reacID <- model$reacID[bs]
        newmodel$namesSpecies <- model$namesSpecies
      }
    } else {
      newmodel <- model
    }
    return(newmodel)
  }


  bin2dec <- function(x) {
    exp2 <- 2^c((length(x)-1):0)
    y <- exp2%*%x
    return(y)
  }

  dec2bin <- function(x) {
    if (x == 0) {
      y <- 0
    } else {
      tmp <- rev(as.integer(intToBits(x)))
      y <- tmp[which(tmp != 0)[1]:length(tmp)] # this is faster than my dec2bin
    }
    return(y)
  }

  dec2binOld <- function(x) {
    if (x == 0) {
      y <- 0
    } else {
      xTmp <- x
      start <- floor(log(xTmp)/log(2))
      y <- numeric(start)
      count <- 0
      for (i in start:0) {
        count <- count + 1
        exp <- floor(log(xTmp)/log(2))
        if (i == exp) {
          y[count] <- 1
          xTmp <- xTmp - 2^exp
        } else {
          y[count] <- 0
        }
      }
    }
    return(y)
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
      sfInit(parallel=TRUE, cpus=parallel)
    }
    sfLibrary(CellNOptR)
    sfExport(list=exportVars("ex"), local = T)
  }
  spaceExp <- 2^length(model$reacID)
  print(paste("structurally different networks: ", spaceExp, sep = ""))
  ## reduce to equivalent classes by absorption
  bigBang <- function(j) {
    ## cat(".")
    essential <- dec2bin((j-1))
    ## tmp <- rev(as.integer(intToBits(j)))
    ## essential <- tmp[which(tmp != 0)[1]:length(tmp)] # this is faster than my dec2bin
    ## tmp <- absorption(c(rep(0, (length(model$reacID)-length(essential))), essential), model)
    tmp <- reduceGraph(c(rep(0, (length(model$reacID)-length(essential))), essential), model, CNOlist)
    return(tmp)
  }
  pop <- matrix(NA, nrow = spaceExp, ncol = length(model$reacID))
  if (!is.null(parallel)) {
    pop <- sfLapply(as.list(1:spaceExp), bigBang)
  } else {
    pop <- lapply(as.list(1:spaceExp), bigBang)
    cat("\n")
  }
  pop <- do.call("rbind", pop)
  res <- pop
  res <- apply(res, 1, paste, collapse = "")
  realSpace <- 1:spaceExp
  if (sum(duplicated(res) == TRUE) > 0 & reduce) {
    realSpace <- realSpace[-which(duplicated(res) == TRUE)]
  }
  ## pop <- pop[-which(duplicated(res) == TRUE), ]
  print(paste("reduction to structurally different equivalence classes: ", length(realSpace), sep = ""))
  scores <- numeric(spaceExp)
  getBinScore <- function(j) {
    cat(".")
    ## essential <- dec2bin((j-1))
    ## pop[j, ] <- c(rep(0, (ncol(pop)-length(essential))), essential)
    scores[j] <- computeScoreNemT1(CNOlist, model = model, pop[j, ], sizeFac = sizeFac, NAFac = NAFac, NEMlist = NEMlist, parameters = parameters, method = method, relFit = relFit)
    return(scores[j])
  }
  if (!is.null(parallel)) {
    res <- sfApply(as.matrix(realSpace), 1, getBinScore)
  } else {
    res <- apply(as.matrix(realSpace), 1, getBinScore)
    cat("\n")
  }
  ## samples <- character()
  ## for (j in realSpace) {
  ##   essential <- dec2bin((j-1))
  ##   pop[j, ] <- c(rep(0, (ncol(pop)-length(essential))), essential)
  ##   samples <- c(samples, toString(pop[j, ]))
  ## }
  bString <- (pop[realSpace, ])[which(res == min(res))[1], ]
  dtmRatio <- mean(abs(res - mean(res)))/abs(min(res) - mean(res))
  dtmRatio2 <- 1 - abs(min(res) - mean(res))/abs(min(res) - max(res))
  print("best network found:")
  print(toString(bString))
  print("score:")
  print(min(res))
  print("distance to mean ratio (smaller is better):")
  print(dtmRatio)
  print(dtmRatio2)
  #hist(res)
  if (!is.null(parallel)) {
    sfStop()
  }
  plotDnf(model$reacID[as.logical(bString)])
  ## names(res) <- samples
  return(list(bString = bString, score = min(res), bStrings = pop, scores = res, dtmRatio = dtmRatio))
}
