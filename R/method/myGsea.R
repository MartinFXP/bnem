myGsea <- function(testList, goList, parallel = NULL, adjust.method = "FDR", conservative = TRUE) {
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
    sfExport("testList", "goList")
  }
  
  startTime <- Sys.time()

  if (conservative) {
    Complete <- unique(intersect(unlist(testList),unlist(goList))) # use this for conservative method
  } else {
    Complete <- unique(c(unlist(testList),unlist(goList))) # use this for none conservative method
  }
  N <- length(Complete)
  
  checkGeneCluster <- function(j,i) {

    resDf <- data.frame()
    genelist <- unique(intersect(testList[[i]],Complete))
    
    if (!is.null(genelist)) {
      targetlist <- intersect(goList[[j]],Complete)
      n <- length(targetlist)
      D <- length(genelist)
      AinB <- length(intersect(targetlist, genelist))
      fisherMat <- cbind(c(AinB,D-AinB),c(n-AinB,N-n-D+AinB))
      fisherTest <- fisher.test(fisherMat, alternative = "greater")
      ABpval <- fisherTest$p.value
      resDf <- rbind(resDf, data.frame(test = names(testList)[i], go = names(goList)[j], pval = ABpval, overlap = AinB, targetset = n, testset = D))  
    }
    
    return(resDf)
  }
  
  dfList <- NULL
  if (is.null(parallel)) {
    for (i in 1:length(testList)) {
      cat(".")
      dfList <- c(dfList, lapply(1:length(goList), checkGeneCluster, i))
    }
  } else { 
    for (i in 1:length(testList)) {
      dfList <- c(dfList, sfLapply(1:length(goList), checkGeneCluster, i))
    }
  }

  cat("\n")
  print(Sys.time() - startTime)
  print("all hypergeometric tests (=fisher tests with 'greater') completed...")
  print("preparing data...")
  
  bigDf <- do.call("rbind", dfList)
  
  bigDf <- bigDf[order(bigDf[, 3]), ]

  if (adjust.method %in% "FDR") {

    qvals <- (nrow(bigDf)*bigDf[, 3])/(1:nrow(bigDf))

    for (i in (length(qvals)-1):1) {

      qvals[i] <- min(qvals[i], qvals[i+1])

    }

    qvals[which(qvals > 1)] <- 1 # except for this line identical to benjamini-hochberg

  } else {
    
    qvals <- p.adjust(bigDf[, 3], method = adjust.method)

  }
  
  bigDf <- cbind(bigDf, qvals = qvals)

  bigDf <- bigDf[order(bigDf[, 7]), ]
  
  if (!is.null(parallel)) {
    sfStop()
  }

  return(bigDf)
}
