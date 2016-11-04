resBNEM <-
function(bString, CNOlist, model, NEMlist, parameters = list(cutOffs = c(0,1,0), scoring = c(0.1,0.2,0.9)), method = "s", sizeFac = 10^-10, main = "residuals for decoupled vertices", sub = "green residuals are added effects (left positive, right negative) and red residuals are deleted effects", cut = TRUE, approach = "fc", parallel = NULL, verbose = TRUE, ...) {
  CNOlist <- checkCNOlist(CNOlist)
  method <- checkMethod(method)[1]
  NEMlist <- checkNEMlist(NEMlist, CNOlist, parameters, approach, method)
  simResults <- simulateStatesRecursive(CNOlist = CNOlist, model = model, bString = bString)
  SCompMat <- computeFc(CNOlist, t(simResults))
  SCompMat <- SCompMat[, which(colnames(SCompMat) %in% colnames(NEMlist$fc))]
  NEMlist$fc <- NEMlist$fc[, order(colnames(NEMlist$fc))]
  SCompMat <- SCompMat[, order(colnames(SCompMat))]
  SCompMat <- SCompMat[, colnames(NEMlist$fc)]
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- colnames(CNOlist@inhibitors)
  tmp <- computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 1, parameters = parameters, method = method, verbose = verbose, sizeFac = sizeFac)
  EtoS <- tmp$EtoS

  if (verbose) {
    print(paste("calculating residuals for ", ncol(CNOlist@signals[[1]]), " S-genes based on ", length(unique(EtoS[, 1])), " E-genes.", sep = ""))
  }
  
  resMat <- matrix(0, nrow = ncol(CNOlist@signals[[1]]), ncol = 2*ncol(NEMlist$fc))
  resVec <- numeric(ncol(CNOlist@signals[[1]]))
  resType <- matrix(0, nrow = ncol(CNOlist@signals[[1]]), ncol = 2*ncol(NEMlist$fc))

  checkSgene <- function(i) {
    resType <- numeric(2*ncol(NEMlist$fc))
    resMat <- numeric(2*ncol(NEMlist$fc))
    if (sum(EtoS[, 2] == i) == 0) {
      resType[1:length(resType)] <- -1
      resMat[1:length(resType)] <- -1
    } else {
      if (sum(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))) == 1) {
        data.tmp <- t(as.matrix(NEMlist$fc[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))), ]))
        rownames(data.tmp) <- rownames(NEMlist$fc)[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i)))]
      } else {
        data.tmp <- NEMlist$fc[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))), ]
      }
      resVec <- sum(abs(cor(SCompMat[i, ], t(data.tmp), method = method)))
      for (j in 1:ncol(data.tmp)) { # parallel this!
        
        cat('\r', paste(floor(((i-1)*ncol(data.tmp) + j)/(ncol(CNOlist@signals[[1]])*ncol(data.tmp))*100), "%", sep = ""))
        flush.console()
        
        sgene <- SCompMat[i, ]
        mem <- sgene[j]
        if (mem == 0) {
          resType[j] <- 1
          resType[(j+ncol(data.tmp))] <- 1
          sgene[j] <- 1
          resMat[j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
          sgene[j] <- -1
          resMat[(j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
        }
        if (mem == 1) {
          resType[j] <- -1
          resType[(j+ncol(data.tmp))] <- 1
          sgene[j] <- 0
          resMat[j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
          sgene[j] <- -1
          resMat[(j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
        }
        if (mem == -1) {
          resType[j] <- 1
          resType[(j+ncol(data.tmp))] <- -1
          sgene[j] <- 1
          resMat[j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
          sgene[j] <- 0
          resMat[(j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
        }
      }
    }
    return(list(resMat = resMat, resType = resType, resVec = resVec))
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
      sfInit(parallel=TRUE, cpus=parallel, type = "SOCK")
    }
    sfLibrary(CellNOptR)
    library(bnem)
    ## sfExport(list = c("checkSgene", "computeScoreNemT1", "simulateStatesRecursiveAdd", "simulateStatesRecursive", "reduceGraph", "getNemFit", "checkCNOlist", "checkNEMlist", "computeFc",  "computeSm", "sizeFac", "method", "removeCycles", "dnf2adj", "plotBinary", "getHierarchy", "absorption", "checkMethod", "approach", "parameters"), local = T)
  }
  
  if (!is.null(parallel)) {
    resTmp <- sfLapply(as.list(1:nrow(resMat)), checkSgene)
    sfStop()
  } else {
    resTmp <- lapply(as.list(1:nrow(resMat)), checkSgene)
  }
  
  for (i in 1:nrow(resMat)) {
    resMat[i, ] <- resTmp[[i]]$resMat
    resType[i, ] <- resTmp[[i]]$resType
    resVec[i] <- resTmp[[i]]$resVec
  }
  
  ## resMat <- do.call("rbind", resTmp)
    
  ## for (i in 1:ncol(CNOlist@signals[[1]])) {
  ##   if (sum(EtoS[, 2] == i) == 0) {
  ##     resType[i, ] <- -1
  ##     resType[i, ] <- -1
  ##     next()
  ##   }
  ##   if (sum(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))) == 1) {
  ##     data.tmp <- t(as.matrix(NEMlist$fc[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))), ]))
  ##     rownames(data.tmp) <- rownames(NEMlist$fc)[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i)))]
  ##   } else {
  ##     data.tmp <- NEMlist$fc[which(rownames(NEMlist$fc) %in% names(which(EtoS[, 2] == i))), ]
  ##   }
  ##   resVec[i] <- sum(abs(cor(SCompMat[i, ], t(data.tmp), method = method)))
  ##   for (j in 1:ncol(data.tmp)) {
      
  ##     cat('\r', paste(floor(((i-1)*ncol(data.tmp) + j)/(ncol(CNOlist@signals[[1]])*ncol(data.tmp))*100), "%", sep = ""))
  ##     flush.console()
      
  ##     sgene <- SCompMat[i, ]
  ##     mem <- sgene[j]
  ##     if (mem == 0) {
  ##       resType[i, j] <- 1
  ##       resType[i, (j+ncol(data.tmp))] <- 1
  ##       sgene[j] <- 1
  ##       resMat[i, j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##       sgene[j] <- -1
  ##       resMat[i, (j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##     }
  ##     if (mem == 1) {
  ##       resType[i, j] <- -1
  ##       resType[i, (j+ncol(data.tmp))] <- 1
  ##       sgene[j] <- 0
  ##       resMat[i, j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##       sgene[j] <- -1
  ##       resMat[i, (j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##     }
  ##     if (mem == -1) {
  ##       resType[i, j] <- 1
  ##       resType[i, (j+ncol(data.tmp))] <- -1
  ##       sgene[j] <- 1
  ##       resMat[i, j] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##       sgene[j] <- 0
  ##       resMat[i, (j+ncol(data.tmp))] <- sum(abs(cor(sgene, t(data.tmp), method = method)))
  ##     }
  ##   }
  ## }
  resType <- resType*(-1)
  resDiff <- (resVec - resMat)/nrow(NEMlist$fc)
  colnames(resDiff) <- c(colnames(NEMlist$fc), colnames(NEMlist$fc))
  rownames(resDiff) <- colnames(CNOlist@signals[[1]])
  resDiff1 <- cbind(resDiff[, 1:ncol(NEMlist$fc)], max(resDiff), resDiff[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])

  p1 <- heatmapOP(resDiff1, Rowv = F, Colv = F, main = main, sub = sub, bordercol = "grey", ...)

  resDiff2 <- cbind(resDiff[, 1:ncol(NEMlist$fc)], min(resDiff), resDiff[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])
  resType2 <- cbind(resType[, 1:ncol(NEMlist$fc)], min(resType), resType[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])
  resDiff2[which(resDiff2 > 0)] <- 0

  p2 <- heatmapOP(resDiff2, Colv = F, Rowv = F, main = main, sub = sub, bordercol = "grey", ...)

  resDiff3 <- resDiff2*resType2

  p3 <- heatmapOP(resDiff3, Colv = F, Rowv = F, main = main, sub = sub, bordercol = "grey", ...)

  #print(p1, position=c(0, 0, .33, 1), more=TRUE)
  #print(p2, position=c(.33, 0, .66, 1), more=TRUE)
  #print(p3, position=c(.66, 0, 1, 1))

  res.breaks <- seq(-max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))), max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))), (max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))) - -max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))))/100)

  p1 <- heatmapOP(resDiff3[, 1:ncol(NEMlist$fc)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (positive effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = FALSE)
  
  p2 <- heatmapOP(resDiff3[, (ncol(NEMlist$fc)+2):(2*ncol(NEMlist$fc)+1)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (negative effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = T)

  if (verbose) {
    print(p1, position=c(0, 0, .48, 1), more=TRUE)
    print(p2, position=c(.48, 0, 1, 1))
  }
  
  #print(p3)

  if (cut & all(is.na(resDiff) == F)) {
    if (sum(apply(abs(resDiff1), 1, sum) == 0) > 0) {
      resDiff1 <- resDiff1[-which(apply(abs(resDiff1), 1, sum) == 0), ]
    }
    if (sum(apply(abs(resDiff1), 2, sum) == 0) > 0) {
      resDiff1 <- resDiff1[, -which(apply(abs(resDiff1), 2, sum) == 0)]
    }
    if (sum(apply(abs(resDiff2), 1, sum) == 0) > 0) {
      resDiff2 <- resDiff2[-which(apply(abs(resDiff2), 1, sum) == 0), ]
    }
    if (sum(apply(abs(resDiff2), 2, sum) == 0) > 0) {
      resDiff2 <- resDiff2[, -which(apply(abs(resDiff2), 2, sum) == 0)]
    }
    if (sum(apply(abs(resDiff3), 1, sum) == 0) > 0) {
      resDiff3 <- resDiff3[-which(apply(abs(resDiff3), 1, sum) == 0), ]
    }
    if (sum(apply(abs(resDiff3), 2, sum) == 0) > 0) {
      resDiff3 <- resDiff3[, -which(apply(abs(resDiff3), 2, sum) == 0)]
    }
  }
  return(list(resDiff1 = resDiff1, resDiff2 = resDiff2, resDiff3 = resDiff3))
}
