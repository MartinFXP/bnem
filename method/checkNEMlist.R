checkNEMlist <- function(NEMlist, CNOlist, parameters, approach, method) {
  NEMlistTmp <- NEMlist
  if("abs" %in% approach) {
    if (length(table(NEMlist$exprs)) == 2) {
      NEMlist$norm <- NEMlist$exprs
      NEMlist$norm[which(NEMlist$norm == as.numeric(names(table(NEMlist$exprs))[1]))] <- 0
      NEMlist$norm[which(NEMlist$norm == as.numeric(names(table(NEMlist$exprs))[2]))] <- 1
    }
    if (length(NEMlist$norm) == 0 & !any(c("pearson", "spearman", "kendall") %in% method)) {
      if ("pam" %in% approach) {
        print("data is not discretized/normed to (0,1); performing two cluster pam normalization")
        NEMlist$norm <- pamNorm(NEMlist$exprs)
      } else {
        print("data is not discretized/normed to (0,1); performing simple normalization")
        NEMlist$norm <- simpleNorm(NEMlist$exprs)
      }
      if ("kmeans" %in% approach) {
        print("data is not discretized/normed to (0,1); performing two means normalization")
        NEMlist$norm <- pamNorm(NEMlist$exprs)
      } else {
        print("data is not discretized/normed to (0,1); performing simple normalization")
        NEMlist$norm <- simpleNorm(NEMlist$exprs)
      }
    }
    if ("spearman" %in% method) {
      colnames.exprs <- colnames(NEMlist$exprs)
      rownames.exprs <- rownames(NEMlist$exprs)
      NEMlist$exprs <- t(apply(NEMlist$exprs, 1, rank))
      colnames(NEMlist$exprs) <- colnames.exprs
      rownames(NEMlist$exprs) <- rownames.exprs
    }
  }
  if ("fc" %in% approach) {
    if (length(NEMlist$fc) == 0)  {
      print("foldchanges missing; automatic calculation")
      NEMlist$fc <- computeFc(CNOlist, NEMlist$exprs) 
      egenes <- NEMlist$egenes
      if (parameters$cutOffs[2] != 0) {
        NEMlist <- computeSm(CNOlist, NEMlist, parameters, method = method)
      }
      NEMlist$egenes <- egenes
    }
    if (length(NEMlist$E0) == 0)  {
      egenes <- NEMlist$egenes
      if (parameters$cutOffs[2] != 0) {
        NEMlist <- computeSm(CNOlist, NEMlist, parameters, method = method)
      }
      NEMlist$egenes <- egenes
    }
    if ("spearman" %in% method) {
      colnames.fc <- colnames(NEMlist$fc)
      rownames.fc <- rownames(NEMlist$fc)
      NEMlist$fc <- t(apply(NEMlist$fc, 1, rank))
      colnames(NEMlist$fc) <- colnames.fc
      rownames(NEMlist$fc) <- rownames.fc
    }
  }
  if (!is.null(NEMlist$egenes) & is.null(NEMlist$geneGrid)) {
    sgeneAdd <- matrix(Inf, nrow = nrow(NEMlist$fc), ncol = (ncol(CNOlist@signals[[1]])*2))
    for (i in 1:nrow(NEMlist$fc)) {
      egeneName <- rownames(NEMlist$fc)[i]
      sgeneCols <- numeric()
      for (j in 1:length(NEMlist$egenes)) {
        if (egeneName %in% NEMlist$egenes[[j]]) {
          colTmp <- which(colnames(CNOlist@signals[[1]]) == names(NEMlist$egenes)[j])
          sgeneCols <- c(sgeneCols, colTmp, (colTmp+ncol(CNOlist@signals[[1]])))
        }
      }
      sgeneAdd[i, sgeneCols] <- 0
    }
    sgeneAddCheck <- rowMins(sgeneAdd)
    sgeneAdd[which(sgeneAddCheck == Inf), ] <- 0
    NEMlist$geneGrid <- sgeneAdd
  }
  if (!is.null(NEMlistTmp$weights)) {
    NEMlist$weights <- NEMlistTmp$weights
  }
  NEMlist$signalStates <- NEMlistTmp$signalStates
  return(NEMlist)
}
