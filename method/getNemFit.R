getNemFit <- function (simResults, CNOlist, model, indexList, # VERSION of CNO: 1.4
                       timePoint = c("t1", "t2"),
                       NAFac = 1, sizePen, simResultsT0 = NA,
                       tellme = 0,
                       approach = "fc",
                       NEMlist,
                       parameters,
                       sim = 0,
                       relFit = FALSE,
                       method = "pearson",
                       verbose = FALSE,
                       opt = "min"
                       ) {
  if ((class(CNOlist) == "CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  ## if ((length(bString) <= 100 & sim == 0) | sim == 2) {
  ##   simResults <- simResults[, indexList$signals] # not necessary with new simulatestates function
  ## }
  ## if (length(bString) > 100 | sim == 1) {
  #simResults <- simResults[, which(colnames(simResults) %in% colnames(CNOlist@signals[[1]]))] # but that is
  simResults <- simResults[, colnames(CNOlist@signals[[1]])] # but that is
  ## }
  ## keep that stuff if I want to look at different timepoints:
  if (timePoint == "t1") {
    tPt <- 2
  }
  else {
    if (timePoint == "t2") {
      tPt <- 3
    }
    else {
      tPt <- timePoint
    }
  }
  # MY CODE START start with the nem approaches for the error
  MSEAabs <- NULL
  MSEIabs <- NULL
  MSEAfc <- NULL
  MSEIfc <- NULL
  require(matrixStats)
  CNOlist <- checkCNOlist(CNOlist)
  NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist, parameters = parameters, approach = approach, method = method)
  if ("abs" %in% approach) {
    MSEE <- rep(Inf, nrow(NEMlist$norm)) # the mse for each egene for the current simulated steady state given the current hyperedges
    S.mat <- simResults
    colnames(NEMlist$norm)[which(colnames(NEMlist$norm) %in% "")] <- "Base"
    rownames(S.mat)[which(rownames(S.mat) %in% "")] <- "Base"
    S.mat <- S.mat[colnames(NEMlist$norm), ]
    if (any(c("spearman", "pearson", "kendall") %in% method)) {
      if ("spearman" %in% method) {
        E.mat <- NEMlist$exprs
        S.mat.ranks <- t(apply(S.mat, 1, rank))
        cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p"))
      } else {
        E.mat <- NEMlist$exprs
        cosine.sim <- -t(cor(S.mat, t(E.mat), method = method))
      }
      R <- cbind(cosine.sim, -cosine.sim)
    } else {
      ESpos <- NEMlist$norm%*%abs(1 - S.mat)*parameters$scoring[1] # count mismatches for zeros
      ES0 <- abs(1 - NEMlist$norm)%*%S.mat*parameters$scoring[2] # count mismatches for ones
      ESposI <- NEMlist$norm%*%S.mat*parameters$scoring[1] # count mismatches for zeros if inhibited
      ES0I <- abs(1 - NEMlist$norm)%*%abs(1 - S.mat)*parameters$scoring[2] # count mismatches for ones if inhibited
      MSEAabs <- (ESpos + ES0)
      MSEIabs <- (ESposI + ES0I)
      R <- cbind(MSEAabs, MSEIabs)
      R <- R/ncol(NEMlist$exprs)
    }
    R[is.na(R)] <- max(R[!is.na(R)])
    MSEE <- rowMins(R)
  }
  if ("fc" %in% approach) {
    MSEE <- rep(Inf, nrow(NEMlist$fc))
    SCompMat <- computeFc(CNOlist, t(simResults))
    SCompMat[SCompMat > 0] <- 1
    SCompMat[SCompMat < 0] <- -1
    SCompMat <- t(SCompMat)
    ## print(colnames(NEMlist$fc)[which(!(colnames(NEMlist$fc) %in% rownames(SCompMat)))]); print(rownames(SCompMat)[which(!(rownames(SCompMat) %in% colnames(NEMlist$fc)))]) # for debugging
    ## make Sgene matrix and Egene matrix foldchanges have the same number and/or order of experiments:
    if (is.null(rownames(SCompMat))) {
      SCompMat <- SCompMat[, colnames(NEMlist$fc)]
    } else {
      SCompMat <- SCompMat[colnames(NEMlist$fc), ] # replicates now possible? that was easy!
    }
    ## NOTE: either check the 0s also or think of some standardweight to penalize networks where only a few changes occur
    if (any(c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "spearman", "pearson", "kendall") %in% method)) {
      if (ncol(SCompMat) != ncol(NEMlist$fc)) {
        SCompMat <- t(SCompMat)
      }
      if (length(grep("_vs_", rownames(SCompMat))) > 0) {
        SCompMat <- t(SCompMat)
      }
      S.mat <- SCompMat
      E.mat <- NEMlist$fc
      if (any(c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") %in% method)) {
        power <- as.numeric(method)
        power <- power[-which(is.na(power)==T)]
        require(flexclust)
        if (length(power) == 0) {
          power <- 2
        }
        method <- method[which(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))[1]]
        cosine.sim <-  dist2(S.mat, E.mat, method = method, p = power)
        cosine.simI <- dist2(S.mat, -E.mat, method = method, p = power)
        MSEAfc <- t(cosine.sim)
        MSEIfc <- t(cosine.simI)
      }
      if (any(c("spearman", "pearson", "kendall", "cosine") %in% method)) {
        weighted <- parameters$scoring[2]
        if (weighted > 1) {
          for (i in 1:nrow(S.mat)) {
            S.mat[i, grep(rownames(S.mat)[i], colnames(S.mat))] <- S.mat[i, grep(rownames(S.mat)[i], colnames(S.mat))]*weighted
          }
        }
        require(matrixStats)
        ## if (all(S.mat==0)) {
        ##   noise <- rnorm(ncol(S.mat), 0, 0.00001) # that is a really stupid idea
        ##   S.mat <- t(t(S.mat) + noise)
        ##   E.mat <- t(t(E.mat) + noise)
        ## }
        if (parameters$scoring[1] == 0) {
          S.mat[which(S.mat == 0)] <- NA
        }
        if ("pearson" %in% method) {
          cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p", use = "pairwise.complete.obs"))
        }
        if ("spearman" %in% method) {
          S.mat.ranks <- apply(S.mat, 1, rank)
          cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p", use = "pairwise.complete.obs"))
        }
        if ("kendall" %in% method) {
          cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "k", use = "pairwise.complete.obs"))
        }
        if ("cosine" %in% method) {
          if (any(is.na(S.mat))) {
            S.mat[is.na(S.mat)] <- 0
          }
          vprod <- function(x) { return(x%*%x) }
          S.mat <- S.mat/(apply(S.mat, 1, vprod)^0.5)
          E.mat <- E.mat/(apply(E.mat, 1, vprod)^0.5)
          cosine.sim <- -t(S.mat%*%t(E.mat))
        }
        if ("test" %in% method) {
          cs.fisher <- atanh(cosine.sim) # 0.5*log((1+cosine.sim)/(1-cosine.sim))
          cs.sd <- 1/sqrt(ncol(NEMlist$fc) - 3)
          cs.z <- (cs.fisher - 0)/cs.sd
          cs.sign <- sign(cosine.sim)
          cosine.sim <- 2*pnorm(-abs(cs.z)) # low p-value means high score
          cosine.sim[1:length(cosine.sim)] <- abs(1 - p.adjust(cosine.sim, method = "bonferroni"))
          cosine.sim <- cosine.sim*cs.sign
        }
        cosine.sim[is.na(cosine.sim)] <- 0
        MSEAfc <- cosine.sim
        MSEIfc <- -cosine.sim
      }
    } else {
      S0 <- as.matrix(1 - abs(SCompMat))
      Spos <- SCompMat
      Spos[which(Spos == -1)] <- 0
      Sneg <- SCompMat
      Sneg[which(Sneg == 1)] <- 0
      if ("mLL" %in% method | "cp" %in% method) {
        Sneg <- abs(Sneg)
        E0 <- NEMlist$E0
        Epos <- NEMlist$Epos
        Eneg <- NEMlist$Eneg
        EposI <- NEMlist$EposI
        EnegI <- NEMlist$EnegI
        ES0 <- E0%*%S0
        ES0pos <- E0%*%Spos
        ES0neg <- E0%*%Sneg
        ESpos <- Epos%*%Spos
        ESpos0 <- Epos%*%S0
        ESposneg <- Epos%*%Sneg
        ESneg <- Eneg%*%Sneg
        ESneg0 <- Eneg%*%S0
        ESnegpos <- Eneg%*%Spos
        ESposI <- EposI%*%Spos
        ESposI0 <- EposI%*%S0
        ESposIneg <- EposI%*%Sneg
        ESnegI <- EnegI%*%Sneg
        ESnegI0 <- EnegI%*%S0
        ESnegIpos <- EnegI%*%Spos
        alpha1 <- parameters$scoring[1]
        beta2 <- parameters$scoring[2]
        beta4 <- beta2/3
        if ("cont" %in% method) {
          MSEAfc <- exp(ESpos + ESneg)
          MSEIfc <- exp(ESposneg + ESnegpos)
        } else {
          alpha2 <- alpha1/2
          beta1 <- parameters$scoring[3]
          beta3 <- beta2/2
          Alpha <- 1 - alpha1 - alpha2
          Beta <- 1 - beta1 - beta2 - beta3 - beta4
          MSEAfc <- Alpha^ES0 * Beta^(ESpos + ESneg) * beta1^(ESposI + ESnegI) * beta2^(ES0pos + ES0neg) * beta3^(ESposIneg + ESnegIpos) * beta4^(ESposneg + ESnegpos) * alpha1^(ESposI0 + ESnegI0) * alpha2^(ESpos0 + ESneg0)
          MSEIfc <- Alpha^ES0 * beta4^(ESpos + ESneg) * beta3^(ESposI + ESnegI) * beta2^(ES0pos + ES0neg) * beta1^(ESposIneg + ESnegIpos) * Beta^(ESposneg + ESnegpos) * alpha1^(ESposI0 + ESnegI0) * alpha2^(ESpos0 + ESneg0)
        }
        if ("cp" %in% method) {
          MSEAfc <- -log(MSEAfc)
          MSEIfc <- -log(MSEIfc)
        }
      } else {
        E0 <- NEMlist$E0
        Epos <- NEMlist$Epos
        Eneg <- NEMlist$Eneg
        EposI <- NEMlist$EposI
        EnegI <- NEMlist$EnegI
        ## differ between count all fc and alignment
        ESpos <- Epos%*%Spos
        ESneg <- Eneg%*%Sneg
        ESposI <- EposI%*%Spos
        ESnegI <- EnegI%*%Sneg
        ES0 <- E0%*%S0
        MSEAfc <- (ESpos + ESneg + ES0)
        MSEIfc <- (ESposI + ESnegI + ES0)
      }
    }
    if ("mLL" %in% method) {
        R <- rbind(MSEAfc, MSEIfc)
        MSEE <- log(rowSums(R))
    } else {
      if (any(c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "spearman", "pearson", "kendall") %in% method)) {
        R <- cbind(MSEAfc, MSEIfc)
        R[is.na(R)] <- max(R[!is.na(R)])
        if (relFit) {
          S0 <- as.matrix(1 - abs(SCompMat))
          Spos <- SCompMat
          Spos[which(Spos == -1)] <- 0
          Sneg <- SCompMat
          Sneg[which(Sneg == 1)] <- 0
          SS0 <- t(S0*parameters$scoring[1])%*%S0
          SSpos <- t(Spos)%*%Spos
          SSneg <- t(Sneg)%*%Sneg
          MSES <- (SSpos + SSneg + SS0)/ncol(NEMlist$fc)
          R <- t(t(R)*diag(MSES))
        }
      } else {
        R <- cbind(MSEAfc, MSEIfc)
        R[is.na(R)] <- max(R[!is.na(R)])
        if (relFit) {
          SS0 <- t(S0*parameters$scoring[1])%*%S0
          SSpos <- t(Spos)%*%Spos
          SSneg <- t(Sneg)%*%Sneg
          MSES <- (SSpos + SSneg + SS0)
          R <- t(t(R)/diag(MSES))
          R[is.na(R)] <- 0
        } else {
          R <- R/ncol(NEMlist$fc)
        }
      }
      if (is.null(NEMlist$egenes)) {
        MSEE <- rowMins(R)
      } else {
        MSEE <- numeric()
        R <- R+NEMlist$geneGrid
        R[is.na(R)] <- Inf
        if (is.null(NEMlist$weights)) {
          MSEE <- rowMins(R)
        } else {
          R[is.infinite(R)] <- 0
          topNsum <- function(x, N) {
            y <- sum(x[order(x)[1:N]])
            return(y)
          }
          MSEE <- apply(R, 1, topNsum, ncol(CNOlist@signals[[1]]))
          MSEE <- MSEE*NEMlist$weights
        }
      }
    }
  }
  if (parameters$cutOffs[3] == -1) {
    require(stats)
    ## do median polish over gene clusters
    data.med <- NEMlist$fc[1:ncol(CNOlist@signals[[2]]), ]*0
    Epos <- which(R == MSEE, arr.ind = T)
    for (i in 1:ncol(CNOlist@signals[[1]])) {
      tmp <- medpolish(rbind(NEMlist$fc[Epos[which(Epos[, 2] == i), 1], ], -NEMlist$fc[Epos[which(Epos[, 2] == (i+ncol(CNOlist@signals[[2]]))), 1], ]), trace.iter=F)
      data.med[i, ] <- tmp$col
    }
    E.mat <- data.med
    E.mat[is.na(E.mat)] <- 0
    tmp <- which(apply(E.mat, 1, sum) != 0)
    E.mat <- as.matrix(E.mat[which(apply(E.mat, 1, sum) != 0), ])
    rownames(E.mat) <- rownames(S.mat)[tmp]
    NEMlist$fc <- E.mat
    S.mat <- SCompMat
    if (parameters$scoring[1] == 0) {
      S.mat[which(S.mat == 0)] <- NA
    }
    if ("pearson" %in% method) {
      cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p", use = "pairwise.complete.obs"))
    }
    if ("spearman" %in% method) {
      S.mat.ranks <- apply(S.mat, 1, rank)
      cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p", use = "pairwise.complete.obs"))
    }
    cosine.sim[is.na(cosine.sim)] <- 0
    MSEAfc <- cosine.sim
    MSEIfc <- -cosine.sim
    R <- cbind(MSEAfc, MSEIfc)
    R[is.na(R)] <- max(R[!is.na(R)])
    MSEE <- rowMins(R)
  }
  if (tellme == 1) {
    names(MSEE) <- rownames(NEMlist$fc)
    if (parameters$cutOffs[3] > 0 & parameters$cutOffs[3] <= 1) {
      R <- R[which(MSEE < -parameters$cutOffs[3]), ]
      MSEE <- MSEE[which(MSEE < -parameters$cutOffs[3])]
    }
    if (parameters$cutOffs[3] < 0 & parameters$cutOffs[3] > -1) {
      score.quantile <- quantile(MSEE, -parameters$cutOffs[3])
      R <- R[which(MSEE < score.quantile), ]
      MSEE <- MSEE[which(MSEE < score.quantile)]
    }
    if (parameters$cutOffs[3] > 1 & parameters$cutOffs[3] <= length(MSEE)) {
      R <- R[order(MSEE)[1:parameters$cutOffs[3]], ]
      MSEE <- MSEE[order(MSEE)[1:parameters$cutOffs[3]]]
    }
    topEgenes <- 0
    subtopo <- matrix(0, nrow = length(MSEE), ncol = ncol(CNOlist@signals[[1]]))
    colnames(subtopo) <- colnames(CNOlist@signals[[1]])
    if (is.null(NEMlist$weights)) {
      Epos <- which(R == MSEE, arr.ind = T)
    } else {
      Epos <- which(NEMlist$geneGrid[, 1:ncol(CNOlist@signals[[1]])] == 0, arr.ind = T)
    }
    if (length(Epos) != 0) {
      if (is.null(dim(Epos))) {
        Epos <- t(as.matrix(Epos))
      }
      posReg <- Epos[which(Epos[, 2] <= ncol(CNOlist@signals[[1]])), ]
      if (length(posReg) > 0) {
        if (length(posReg) > 3) {
          subtopo[posReg] <- 1
        } else {
          subtopo[posReg[1], posReg[2]] <- 1
        }
      }
      negReg <- Epos[which(Epos[, 2] > ncol(CNOlist@signals[[1]])), ]
      if (length(negReg) > 0) {
        if (length(negReg) > 3) {
          negReg[, 2] <- negReg[, 2] - ncol(CNOlist@signals[[1]])
          subtopo[negReg] <- -1
        } else {
          negReg[2] <- negReg[2] - ncol(CNOlist@signals[[1]])
          subtopo[posReg[1], posReg[2]] <- 1
        }
      }

      EtoS <- matrix(0, nrow = nrow(Epos), ncol = 4)
      colnames(EtoS) <- c("Egene", "Sgene", "Type", "MSE")
      rownames(EtoS) <- names(MSEE)[Epos[, 1]]

      Epos[which(Epos[, 2] > ncol(subtopo)), 2] <- Epos[which(Epos[, 2] > ncol(subtopo)), 2] - ncol(subtopo)
      
      EtoS[, 1] <- Epos[, 1]
      EtoS[, 2] <- Epos[, 2]
      EtoS[, 3] <- subtopo[cbind(Epos[, 1], Epos[, 2])]
      EtoS[, 4] <- MSEE[Epos[, 1]]
      
      EtoS <- EtoS[order(EtoS[, 4], decreasing = FALSE), ]
      sgeneScore <- numeric(ncol(subtopo))
    } else {
      EtoS <- matrix(0, nrow = 1, ncol = 4)
      colnames(EtoS) <- c("Egene", "Sgene", "Type", "MSE")
    }
    if (verbose) {
      for (j in 1:ncol(subtopo)) {
        print(paste(j, ".", colnames(simResults)[j], ": ", sum(EtoS[, 2] == j), sep = ""))
        print(paste("Activated: ", sum(EtoS[, 2] == j & EtoS[, 3] == 1), sep = ""))
        print(paste("Inhibited: ", sum(EtoS[, 2] == j & EtoS[, 3] == -1), sep = ""))
        print("Summary Score:")
        print(summary(EtoS[which(EtoS[, 2] == j), 4]))
      }
      dups <- sum(duplicated(rownames(Epos)) == TRUE)
      if (dups > 0) {
        used <- sum(EtoS[-which(duplicated(rownames(Epos)) == TRUE), 2] %in% 1:ncol(subtopo))
      } else {
        used <- nrow(EtoS)
      }
      print(paste("Unique genes used: ", (used), " (", round((used/nrow(NEMlist$fc))*100, 2), " %)", sep = ""))
      print(paste("Duplicated genes: ", dups, sep = ""))
      print("Overall fit:")
      print(summary(EtoS[, 4]))
    }
  }
  if ("mLL" %in% method) {
    deviationPen <- -sum(MSEE)
  } else {
    if (parameters$cutOffs[3] > 1 & parameters$cutOffs[3] <= length(MSEE) & tellme == 0) {
      MSEE <- MSEE[order(MSEE, decreasing = F)[1:parameters$cutOffs[3]]]
    }
    if (parameters$cutOffs[3] > 0 & parameters$cutOffs[3] <= 1 & tellme == 0) {
      MSEE <- MSEE[which(MSEE < -parameters$cutOffs[3])]
    }
    if (parameters$cutOffs[3] < 0 & parameters$cutOffs[3] > -1 & tellme == 0) { # that is just bad and not working for obvious reasons
      score.quantile <- quantile(MSEE, -parameters$cutOffs[3])
      MSEE <- MSEE[which(MSEE < score.quantile)]
    }
    if ("cp" %in% method) {
      deviationPen <- sum(MSEE[!is.na(MSEE)])/nrow(NEMlist$fc) # sum over all the mse for all the egenes
    } else {
      if (parameters$cutOffs[3] > 1 & parameters$cutOffs[3] <= length(MSEE) & tellme == 0) {
        deviationPen <- 1 + sum(MSEE[!is.na(MSEE)])/parameters$cutOffs[3] # sum over all the mse for all the egenes
      } else {
        deviationPen <- 1 + sum(MSEE[!is.na(MSEE)])/nrow(NEMlist$fc) # sum over all the mse for all the egenes
      }
    }
  }
  # END of my code
  NAPen <- NAFac * length(which(is.na(simResults))) # why do I still have that?
  if (opt %in% "min") {
    score <- deviationPen + NAPen + sizePen
  } else {
    score <- 1 - deviationPen - NAPen - sizePen
  }
  if (tellme == 1) {
    return(list(EtoS = EtoS, score = score))
  } else {
    return(score)
  }
}
