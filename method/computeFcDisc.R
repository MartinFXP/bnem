computeFcDisc <- function (x, stimuli = NULL, inhibitors = NULL, batches, runs, method = "edgeR") {
  foldchanges <- numeric()
  foldchangesNames <- character()
  design <- makeDesign(x, stimuli, inhibitors, c(batches[-1], runs[-1]))
  design <- design[, c(stimuli, inhibitors)]
  
  require(edgeR)
  ## # do a generalized linear model
  ## design <- makeDesign(x, stimuli, inhibitors, c(batches[-1], runs[-1]))
  ## design[, 1] <- 1
  ## y <- estimateGLMCommonDisp(x,design)
  ## y <- estimateGLMTrendedDisp(y,design)
  ## y <- estimateGLMTagwiseDisp(y,design)
  ## fit <- glmFit(y,design)
  
  expSums <- apply(design[, c(stimuli, inhibitors)], 1, sum)
  grepCtrls <- which(expSums == 0)
  ctrls <- numeric(nrow(design))
  ctrls[grepCtrls] <- 1
  inhSums <- apply(design[, inhibitors], 1, sum)
  stimSums <- apply(design[, stimuli], 1, sum)
  inhMax <- max(inhSums)
  stimMax <- max(stimSums)
  for (k in 1:inhMax) { # nope use expand.grid ...
    for (i in inhibitors) {
      print(paste("Foldchanges: Ctrl_vs_", i, sep = ""))
      fcTmp <- numeric(nrow(x))
      exps <- design[, i]*2
      exps[which(expSums != k)] <- 0
      if (sum(exps >= 1) == 0) { next() }
      ctrlExps <- ctrls+exps
      xTmp <- x[, which(ctrlExps != 0)]
      ctrlExps <- ctrlExps[which(ctrlExps != 0)]
      y <- DGEList(counts=xTmp,group=ctrlExps)
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      tableTmp <- topTags(et, n = nrow(x)) 
      genesDfUp <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC > 0))
      genesDfUp <- rownames(tableTmp$table)[genesDfUp]
      genesDfDn <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC < 0))
      genesDfDn <- rownames(tableTmp$table)[genesDfDn]
      if (length(genesDfUp) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfUp)] <- 1
      }
      if (length(genesDfDn) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfDn)] <- -1
      }
      foldchanges <- cbind(foldchanges, fcTmp)
      foldchangesNames <- c(foldchangesNames, paste("Ctrl_vs_", i, sep = ""))
    }
  }
  # inhib vs stim+inhib
  for (j in inhibitors) {
    grepInhibs <- intersect(which(expSums == 1), which(design[, j] == 1))
    inhibs <- numeric(nrow(design))
    inhibs[grepInhibs] <- 1
    for (i in stimuli) {
      print(paste("Foldchanges: ", j, "_vs_", j, "_", i, sep = ""))
      fcTmp <- numeric(nrow(x))
      exps <- apply(design[, c(i,j)], 1, sum)
      exps[which(exps == 1)] <- 0
      if (sum(exps >= 1) == 0) { next() }
      inhibExps <- inhibs+exps
      xTmp <- x[, which(inhibExps != 0)]
      inhibExps <- inhibExps[which(inhibExps != 0)]
      y <- DGEList(counts=xTmp,group=inhibExps)
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      tableTmp <- topTags(et, n = nrow(x)) 
      genesDfUp <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC > 0))
      genesDfUp <- rownames(tableTmp$table)[genesDfUp]
      genesDfDn <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC < 0))
      genesDfDn <- rownames(tableTmp$table)[genesDfDn]
      if (length(genesDfUp) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfUp)] <- 1
      }
      if (length(genesDfDn) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfDn)] <- -1
      }
      foldchanges <- cbind(foldchanges, fcTmp)
      foldchangesNames <- c(foldchangesNames, paste(j, "_vs_", j, "_", i, sep = ""))
    }
  }
  # stim vs stim+inhib
  for (j in stimuli) {
    grepStims <- intersect(which(expSums == 1), which(design[, j] == 1))
    stims <- numeric(nrow(design))
    stims[grepStims] <- 1
    for (i in inhibitors) {
      print(paste("Foldchanges: ", j, "_vs_", j, "_", i, sep = ""))
      fcTmp <- numeric(nrow(x))
      exps <- apply(design[, c(i,j)], 1, sum)
      exps[which(exps == 1)] <- 0
      if (sum(exps >= 1) == 0) { next() }
      stimExps <- stims+exps
      xTmp <- x[, which(stimExps != 0)]
      stimExps <- stimExps[which(stimExps != 0)]
      y <- DGEList(counts=xTmp,group=stimExps)
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      tableTmp <- topTags(et, n = nrow(x)) 
      genesDfUp <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC > 0))
      genesDfUp <- rownames(tableTmp$table)[genesDfUp]
      genesDfDn <- intersect(which(tableTmp$table$PValue < 0.05), which(tableTmp$table$logFC < 0))
      genesDfDn <- rownames(tableTmp$table)[genesDfDn]
      if (length(genesDfUp) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfUp)] <- 1
      }
      if (length(genesDfDn) >= 1) {
        fcTmp[which(rownames(x) %in% genesDfDn)] <- -1
      }
      foldchanges <- cbind(foldchanges, fcTmp)
      foldchangesNames <- c(foldchangesNames, paste(j, "_vs_", j, "_", i, sep = ""))
    }
  }
  colnames(foldchanges) <- foldchangesNames
  rownames(foldchanges) <- rownames(x)
  return(foldchanges)
}
