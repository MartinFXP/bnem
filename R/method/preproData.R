preproData <- function(x, # matrix with highdimensional data (samples are columns probes are rows), colnames must be named according to genes (inhibited), batches, stimuli, timepoints and control must denote if there is no inhibition or stimulation
                       controls = c("Stimulicontrol", "Inhibitorcontrol"), # a vector of length; the first entry should be the label of the negative stimulation control and the second entry should be the label for negative inhibition control
                       stimuli = c("Tnf", "Bcr"), # stimulated nodes
                       inhibitors = c("Myc", "Ras"), # inhibited genes
                       signals = "none", # "measured genes": if set to "none", inhibitors above are set as inhibited and measured; if set to node names inhibitors are only set as inhibited and signals are set as measured (where egenes can be connected to)
                       times = "none", # c("0h","1h") timepoints
                       batches = c("none") # names of batches
                       ) {
  x <- as.matrix(x)
  x <- x[, order(colnames(x))]
  ## # throw out low expressed genes:
  ## genes.mean <- apply(x, 1, median)
  ## genes.selection <- rownames(x)[which(genes.mean >= parameters[3])]
  ## x <- x[genes.selection, ]
  ## #controlIndex <- colnames(x)[(1:ncol(x))[-grep(paste(signals, collapse = "|"), colnames(x))]]
  ## # quantile normalization to start:
  ## require(preprocessCore)
  ## y <- normalize.quantiles(x)
  ## colnames(y) <- colnames(x)
  ## rownames(y) <- rownames(x)
  ## x <- y
  # use limma on the data:
  if (batches[1] != "none") {
    require(limma)
    design <- rep(1, ncol(x))
    for (i in batches[-1]) {
      temp <- numeric(ncol(x))
      temp[grep(i, colnames(x))] <- 1 
      design <- cbind(design, temp)
    }
    for (i in stimuli) {
      temp <- numeric(ncol(x))
      temp[grep(i, colnames(x))] <- 1 
      design <- cbind(design, temp)
    }
    for (i in inhibitors) {
      temp <- numeric(ncol(x))
      temp[grep(i, colnames(x))] <- 1 
      design <- cbind(design, temp)
    }
    interactions <- character()
    # only stimuliXgene or also geneXgene interaction???? geneXgene leads to overdefined (or whatever that's called) designmatrix with coefficients not estimable...
    for (i in stimuli) {
      for (j in inhibitors) {
        temp <- numeric(ncol(x))
        temp[intersect(grep(i, colnames(x)), grep(j, colnames(x)))] <- 1
        if (sum(temp == 1) >= 1) {
          interactions <- c(interactions, paste(j, "x", i, sep = ""))
          design <- cbind(design, temp)
        }
      }
    }
    colnames(design) <- c("Ctrl", batches[-1], stimuli, inhibitors, interactions)
    # eliminating time effects like batcheffects is a bit crazy. if the times are not modelled in limma the model will just correct the time effects over all times (mean)
    if (times[1] != "none") { # do the same with times maybe ?
      stop("different times are not implemented yet")
      for (i in times[-1]) {
        temp <- numeric(ncol(x))
        temp[grep(i, colnames(x))] <- 1
        design <- cbind(design, temp)
      }
    }
    colnames(design) <- c("Ctrl", batches[-1], stimuli, inhibitors, interactions, times[-1])
    print("some non estimable coefficients are expected (especially with time series data) and can be ignored; nevertheless all non estimable coefficients should be reviewed by the user")
    fit <- lmFit(x, design)
    fit <- eBayes(fit)
    calcout <- which(colnames(fit$coefficients) %in% batches)
    #calcout <- c(which(colnames(fit$coefficients) %in% batches), which(colnames(fit$coefficients) %in% times))
    coef <- fit$coefficients[, -calcout]
    coef[is.na(coef)] <- 0 # nas shouldn't happen but let's set this to be on the safe side. with interactions and times modelled nas still might/will happen.
    design.sparse <- design[, -calcout]
    xlimma <- t(design.sparse%*%t(coef))
    colnames(xlimma) <- colnames(x)
    xlimma <- xlimma[, grep(batches[1], unique(colnames(xlimma)))] # reduced to the first batch by eliminating batch effects (remember to have each condition available in all batches because otherwise it will be calculated out)
  } else {
    xlimma <- x
    fit <- x
  }
  exprs <- xlimma
  return(list(exprs = exprs))
}
