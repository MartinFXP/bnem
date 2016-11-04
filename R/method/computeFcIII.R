## calculate foldchanges, pvalues and log odds with voom and limma analysis

computeFcIII <- function (x, control, stimuli, inhibitors, batches, runs, sep = "_") {
  NEMlist <- list()
  NEMlist$batches <- batches
  NEMlist$runs <- runs
  NEMlist$sep = sep
  NEMlist$exprs <- x

  require(limma)

  genesMean <- apply(x, 1, median)
  x <- x[order(genesMean, decreasing = TRUE)[1:10000], ]
  x <- x[sort(rownames(x)), ]

  design <- makeDesignFull(x, stimuli, inhibitors, batches, runs)

  design[, 1] <- 0

  design[grep(paste(control[2], ".*", control[1], "|", control[1], ".*", control[2], sep = ""), colnames(x)), 1] <- 1

  xVoom <- voom(x, design, plot = F)
  
  fit <- lmFit(xVoom, design)

  NEMlist$fc <- numeric()
  NEMlist$pvals <- numeric()
  NEMlist$B <- numeric()

  tempNames <- character()
  for (i in c(stimuli, inhibitors)) {
    if (is.na(fit$coefficients[1, i])) { next() }
    contDiff <- paste(i, "-Ctrl", sep = "")
    contrast.matrix <- makeContrasts(contDiff, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    tempNames <- c(tempNames, paste("Ctrl_vs_", i, sep = ""))
    NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(x))$logFC[order(topTable(fit2, n = nrow(x))$ID)])
    NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(x))$adj.P.Val[order(topTable(fit2, n = nrow(x))$ID)])
    NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(x))$B[order(topTable(fit2, n = nrow(x))$ID)])
  }
  colnames(NEMlist$fc) <- c(colnames(NEMlist$fc), tempNames)
  colnames(NEMlist$pvals) <- c(colnames(NEMlist$pvals), tempNames)
  colnames(NEMlist$B) <- c(colnames(NEMlist$B), tempNames)
  rownames(NEMlist$fc) <- topTable(fit2, n = nrow(x))$ID[order(topTable(fit2, n = nrow(x))$ID)]
  rownames(NEMlist$pvals) <- topTable(fit2, n = nrow(x))$ID[order(topTable(fit2, n = nrow(x))$ID)]
  rownames(NEMlist$B) <- topTable(fit2, n = nrow(x))$ID[order(topTable(fit2, n = nrow(x))$ID)]

  for (i in stimuli) {
    for (j in inhibitors) {
      if (is.na(fit$coefficients[1, paste(i, "_", j, sep = "")]) | is.na(fit$coefficients[1, i])) { next() }
      contDiff <- paste(i, "_", j, "-", i, sep = "")
      contrast.matrix <- makeContrasts(contDiff, levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      tempNames <- c(tempNames, paste(i, "_vs_", i, "_", j, sep = ""))
      NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(x))$logFC[order(topTable(fit2, n = nrow(x))$ID)])
      NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(x))$adj.P.Val[order(topTable(fit2, n = nrow(x))$ID)])
      NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(x))$B[order(topTable(fit2, n = nrow(x))$ID)])
      
      contDiff <- paste(i, "_", j, "-", j, sep = "")
      contrast.matrix <- makeContrasts(contDiff, levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      tempNames <- c(tempNames, paste(j, "_vs_", j, "_", i, sep = ""))
      NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(x))$logFC[order(topTable(fit2, n = nrow(x))$ID)])
      NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(x))$adj.P.Val[order(topTable(fit2, n = nrow(x))$ID)])
      NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(x))$B[order(topTable(fit2, n = nrow(x))$ID)])
    }
  }

  colnames(NEMlist$fc) <- tempNames
  colnames(NEMlist$pvals) <- tempNames
  colnames(NEMlist$B) <- tempNames

  NEMlist$fc <- NEMlist$fc[, sort(colnames(NEMlist$fc))]
  NEMlist$pvals <- NEMlist$pvals[, sort(colnames(NEMlist$pvals))]
  NEMlist$B <- NEMlist$B[, sort(colnames(NEMlist$B))]

  return(NEMlist)

}
