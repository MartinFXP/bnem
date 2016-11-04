makeDesignFull <- function(x, stimuli, inhibitors, batches = NULL, runs = NULL, method = "raw") {
  design0 <- makeDesign(x, stimuli, inhibitors, c(batches, runs))
  design <- design0

  ctrlsSum <- apply(design[, -grep(paste(c(runs, batches), collapse = "|"), colnames(design))], 1, sum)
  ctrlsSum <- which(ctrlsSum == 0)
  stimuliDesign <- design[, grep(paste(stimuli, collapse = "|"), colnames(design))]
  inhibitorsDesign <- design[, grep(paste(inhibitors, collapse = "|"), colnames(design))]
  if (length(stimuli) == 1) {
    stimuliDesign <- as.matrix(design[, grep(paste(stimuli, collapse = "|"), colnames(design))])
    colnames(stimuliDesign) <- stimuli
  }
  if (length(inhibitors) == 1) {
    inhibitorsDesign <- as.matrix(design[, grep(paste(inhibitors, collapse = "|"), colnames(design))])
    colnames(inhibitorsDesign) <- inhibitors
  }
  if (is.null(stimuli) == TRUE) {
    stimuliSum <- numeric(nrow(design))
  } else {
    if (length(stimuli) == 1) {
      stimuliSum <- stimuliDesign
    } else {
      stimuliSum <- apply(stimuliDesign, 1, sum)
    }
  }
  if (is.null(inhibitors) == TRUE) {
    inhibitorsSum <- numeric(nrow(design))
  } else {
    if (length(inhibitors) == 1) {
      inhibitorsSum <- inhibitorssDesign
    } else {
      inhibitorsSum <- apply(inhibitorsDesign, 1, sum)
    }
  }
  cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors), collapse = "|"), colnames(design))], 1, sum)
  maxStim <- max(stimuliSum)
  maxKd <- max(inhibitorsSum)
  maxCue <- max(cuesSum)

  grepCtrl <- which(cuesSum == 0)
  grepStims <- intersect(which(stimuliSum != 0), which(inhibitorsSum == 0))
  grepKds <- intersect(which(stimuliSum == 0), which(inhibitorsSum != 0))
  grepStimsKds <- intersect(which(stimuliSum != 0), which(inhibitorsSum != 0))

  design <- numeric()
  designNames <- character()
  design <- rep(0, ncol(x))
  design[grepCtrl] <- 1
  designNames <- "Ctrl"
  
  for (i in grepStims) {
    stimNames <- paste(sort(names(which(stimuliDesign[i, ] >= 1))), collapse = "_")
    if (stimNames %in% designNames) {
      design[i, which(designNames %in% stimNames)] <- 1
    } else {
      design <- cbind(design, rep(0, ncol(x)))
      designNames <- c(designNames, stimNames)
      design[i, which(designNames %in% stimNames)] <- 1
    }
  }

  for (i in grepKds) {
    stimNames <- paste(sort(names(which(inhibitorsDesign[i, ] >= 1))), collapse = "_")
    if (stimNames %in% designNames) {
      design[i, which(designNames %in% stimNames)] <- 1
    } else {
      design <- cbind(design, rep(0, ncol(x)))
      designNames <- c(designNames, stimNames)
      design[i, which(designNames %in% stimNames)] <- 1
    }
  }
  
  for (i in grepStimsKds) {
    stimNames <- paste(c(sort(names(which(inhibitorsDesign[i, ] >= 1))), sort(names(which(stimuliDesign[i, ] >= 1)))), collapse = "_")   # paste(sort(names(which(cbind(stimuliDesign, inhibitorsDesign)[i, ] >= 1))), collapse = "_")
    if (stimNames %in% designNames) {
      design[i, which(designNames %in% stimNames)] <- 1
    } else {
      design <- cbind(design, rep(0, ncol(x)))
      designNames <- c(designNames, stimNames)
      design[i, which(designNames %in% stimNames)] <- 1
    }
  }
  
  if (!is.null(batches))  {
    for (i in batches) {
      if (!is.null(runs)) { 
        for (j in runs) {
          tmp <- numeric(ncol(x))
          tmp[intersect(grep(i, colnames(x)), grep(j, colnames(x)))] <- 1
          if (sum(tmp) != 0) {
            design <- cbind(design, tmp)
            designNames <- c(designNames, paste(sort(c(i, j)), collapse = "_"))
          }
        }
      }
    }
  }
  colnames(design) <- designNames
  if (method %in% "inter") {
    for (i in c(stimuli, inhibitors)) {
      design[, i] <- design0[, i]
    }
  }
  #design <- cbind(design, design0[, grep(batches, colnames(design0))])
  return(design)
}
