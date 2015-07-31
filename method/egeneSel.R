egeneSel <- function(CNOlist, NEMlist,
                     type = "limma", # what criteria to use for selection: limma, silence, sd or c("limma", "sd", ...)
                     parameters = c(log2(1.5), 0.05, 1), # some parameters: 1) the foldchange for "limma" or for the minimal silencing effect for "silence" or the threshold for "variance" (e.g. 2 means a foldchange of 2 in the data has to be there for the egene to be considered) 2) the adjusted p-value needed for "limma" 3) the minimal mean expression for a gene to stay in the dataset (e.g. 2 means if gene E has a lower mean expression than 2 it is thrown out)
                     perGene = 1 # = 1 if you want to know which Egenes belong to which Sgene
                         ) {
  x <- NEMlist$exprs
  # throw out low expressed genes:
  genes.mean <- apply(x, 1, median)
  genes.selection <- which(genes.mean >= parameters[3])
  x <- x[genes.selection, ]
  genes.selection <- list()
  genes.selection$limma <- list()
  genes.selection$silence <- list()
  if ("limma" %in% type) {
    design <- rep(1, ncol(x))
    for (i in colnames(CNOlist@stimuli)) {
      temp <- numeric(ncol(x))
      temp[grep(i, colnames(x))] <- 1 
      design <- cbind(design, temp)
    }
    for (i in colnames(CNOlist@inhibitors)) {
      temp <- numeric(ncol(x))
      temp[grep(i, colnames(x))] <- 1 
      design <- cbind(design, temp)
    }
    ## interactions <- character()
    ## for (i in stimuli) {
    ##   for (j in inhibitors) {
    ##     temp <- numeric(ncol(x))
    ##     temp[intersect(grep(i, colnames(x)), grep(j, colnames(x)))] <- 1
    ##     if (sum(temp == 1) >= 1) {
    ##       interactions <- c(interactions, paste(j, "x", i, sep = ""))
    ##       design <- cbind(design, temp)
    ##     }
    ##   }
    ## }
    colnames(design) <- c("Ctrl", stimuli, inhibitors)
    require(limma)
    fit <- lmFit(x, design)
    fit <- eBayes(fit)
    for (i in colnames(CNOlist@inhibitors)) {
      genes.selection$limma[[i]] <- topTable(fit, coef = i, number = nrow(fit$coefficients))$ID[intersect(which(abs(topTable(fit, coef = i, number = nrow(fit$coefficients))$logFC) >= parameters[1]), which(topTable(fit, coef = i, number = nrow(fit$coefficients))$adj.P.Val <= parameters[2]))]
    }
    ## for (i in stimuli) {
    ##   for (j in inhibitors) {
    ##     temp <- numeric(ncol(x))
    ##     temp[intersect(grep(i, colnames(x)), grep(j, colnames(x)))] <- 1
    ##     if (sum(temp == 1) >= 1) {
    ##       genes.selection <- c(genes.selection, topTable(fit, coef = paste(j, "x", i, sep = ""), number = nrow(fit$coefficients))$ID[intersect(which(abs(topTable(fit, coef = paste(j, "x", i, sep = ""), number = nrow(fit$coefficients))$logFC) >= parameters[1]), which(topTable(fit, coef = paste(j, "x", i, sep = ""), number = nrow(fit$coefficients))$adj.P.Val <= parameters[2]))])
    ##     }
    ##   }
    ## } 
  }
  fc <- NEMlist
  if ("silence" %in% type) {
    fc$fc <- computeFc(CNOlist, NEMlist$exprs)
    fcStimuli <- numeric()
    fcInhibitors <- numeric()
    for (i in colnames(CNOlist@stimuli)) {
      temp <- numeric(ncol(fc$fc))
      temp[grep(i, colnames(fc$fc))] <- 1
      fcStimuli <- cbind(fcStimuli, temp)
    }
    for (i in colnames(CNOlist@inhibitors)) {
      temp <- numeric(ncol(fc$fc))
      temp[grep(i, colnames(fc$fc))] <- 1
      fcInhibitors <- cbind(fcInhibitors, temp)
    }
    maxStimuli <- rowSums(fcStimuli)
    maxInhibitors <- rowSums(fcInhibitors)
    for (i in 1:ncol(CNOlist@stimuli)) {
      genes.selection$silence[[colnames(CNOlist@stimuli)[i]]] <- list()
      stimEffect <- fc$fc[, intersect(which(fcStimuli[, i] == 1), intersect(which(maxStimuli == 1), which(maxInhibitors == 0)))[1]]
      for (j in 1:ncol(CNOlist@inhibitors)) {
        inhibitEffect <- fc$fc[, which(colnames(fc$fc) %in% paste(colnames(CNOlist@stimuli)[i], "vs", colnames(CNOlist@stimuli)[i], colnames(CNOlist@inhibitors)[j], sep = "_"))]
        genes.selection$silence[[colnames(CNOlist@stimuli)[i]]][[colnames(CNOlist@inhibitors)[j]]] <- rownames(fc$fc)[intersect(which(stimEffect > parameters[1]), which(inhibitEffect > parameters[1]))]
        genes.selection$silence[[colnames(CNOlist@stimuli)[i]]][[j]] <- c(genes.selection$silence[[colnames(CNOlist@stimuli)[i]]][[colnames(CNOlist@inhibitors)[j]]], rownames(fc$fc)[intersect(which(stimEffect < -parameters[1]), which(inhibitEffect < -parameters[1]))])
      }
    }
  }
  if ("sd" %in% type) {
    genes.var <- apply(x, 1, sd)
    genes.selection$sd <- rownames(x)[which(genes.var >= parameters[1])]
  }
  if (length(unlist(genes.selection)) == 0) {
    print("no genes selected; selecting all")
    genes.selection <- 1:nrow(xlimma)
  }
  return(genes.selection)
}
