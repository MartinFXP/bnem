dummyCNOlist <- function(stimuli = NULL, inhibitors = NULL, maxStim = 0, maxInhibit = 0, signals = NULL) {
  if (is.null(signals)) {
    signals <- c(stimuli, inhibitors)
  }
  stimn <- length(stimuli)
  inhibn <- length(inhibitors)
  ## do the stim design:
  if (maxStim > stimn) {
    maxStim <- stimn
  }
  if (stimn > 0 & maxStim > 0) {
    mat.size <- 0
    for (i in 1:maxStim) {
      mat.size <- mat.size + choose(stimn, i)
    }
    stimDesign <- matrix(0, mat.size, stimn)
    diag(stimDesign) <- 1
    count <- stimn
    if (maxStim > 1) {
      for (i in 2:maxStim) {
        ## design.tmp <- design[(vertn*(i-2)+1):(vertn*(i-1)), c(2:vertn,1)]
        ## diag(design.tmp) <- 1
        ## design <- rbind(design, design.tmp)
        combs <- t(combn(1:stimn, i))
        for (j in 1:nrow(combs)) {
          count <- count + 1
          stimDesign[count, combs[j, ]] <- 1
        }
      }
    }
    colnames(stimDesign) <- stimuli
  }
  ## do the inhib design:
  inhibn <- length(inhibitors)
  if (maxInhibit > inhibn) {
    maxInhibit <- inhibn
  }
  if (inhibn > 0 & maxInhibit > 0) {
    mat.size <- 0
    for (i in 1:maxInhibit) {
      mat.size <- mat.size + choose(inhibn, i)
    }
    inhibDesign <- matrix(0, mat.size, inhibn)
    diag(inhibDesign) <- 1
    count <- inhibn
    if (maxInhibit > 1) {
      for (i in 2:maxInhibit) {
        ## design.tmp <- design[(vertn*(i-2)+1):(vertn*(i-1)), c(2:vertn,1)]
        ## diag(design.tmp) <- 1
        ## design <- rbind(design, design.tmp)
        combs <- t(combn(1:inhibn, i))
        for (j in 1:nrow(combs)) {
          count <- count + 1
          inhibDesign[count, combs[j, ]] <- 1
        }
      }
    }
    colnames(inhibDesign) <- inhibitors
  }
  ## put them together in combinations:
  if (stimn > 0 & inhibn > 0) {
    if (maxStim > 0 & maxInhibit > 0) {
      design <- matrix(0, nrow(stimDesign)*nrow(inhibDesign), stimn+inhibn)
      for (i in 1:nrow(stimDesign)) {
        design[((i-1)*nrow(inhibDesign) + 1):(i*nrow(inhibDesign)), ] <- cbind(stimDesign[rep(i, nrow(inhibDesign)), ], inhibDesign)
      }
    }
    if (maxStim > 0 & maxInhibit == 0) {
      inhibDesign <- matrix(0, nrow(stimDesign), inhibn)
      design <- cbind(stimDesign, inhibDesign)
    }
    if (maxStim == 0 & maxInhibit > 0) {
      stimDesign <- matrix(0, nrow(inhibDesign), stimn)
      design <- cbind(stimDesign, inhibDesign)
    }
    if (maxStim == 0 & maxInhibit == 0) {
      inhibDesign <- matrix(0, 1, inhibn)
      stimDesign <- matrix(0, 1, stimn)
      design <- cbind(stimDesign, inhibDesign)
    }
    colnames(design) <- c(stimuli, inhibitors)
  }
  colnamesdesign <- colnames(design)
  design <- rbind(cbind(stimDesign, matrix(0, nrow(stimDesign), (ncol(design) - ncol(stimDesign)))), cbind(matrix(0, nrow(inhibDesign), (ncol(design) - ncol(inhibDesign))), inhibDesign), design)
  colnames(design) <- colnamesdesign
  ## make signalmatrix:
  signaln <- length(signals)
  if (signaln > 0) {
    signalData <- matrix(0, nrow(design)+1, signaln)
    colnames(signalData) <- signals
  } else {
    signalData <- matrix(0, nrow(design)+1, 1)
  }
  smult <- nrow(design)/nrow(stimDesign)
  imult <- nrow(design)/nrow(inhibDesign)
  design <- rbind(0, design)
  stimDesign <- design[, 1:ncol(stimDesign)]
  inhibDesign <- design[, (ncol(stimDesign)+1):ncol(design)]
  rownames(design) <- rownames(inhibDesign) <- rownames(stimDesign) <- rownames(signalData) <- c("Ctrl", 2:nrow(design))
  for (i in 2:nrow(design)) {
    tmp <- paste(colnames(design)[which(design[i, ] == 1)], collapse = "_")
    rownames(design)[i] <- rownames(inhibDesign)[i] <- rownames(stimDesign)[i] <- rownames(signalData)[i] <- tmp
  }
  cnolist <- new("CNOlist",
                 cues = design, inhibitors = inhibDesign,
                 stimuli = stimDesign,
                 signals = list(signalData, signalData), timepoints = as.character(c(0,1)))
  cnolist <- checkCNOlist(cnolist)
  return(cnolist)
}

## CNOlist <- dummyCNOlist(paste("S", 1:10, sep = ""), paste("I", 1:10, sep = ""), maxStim = 2, maxInhibit = 1)

## dummyCNOlistOld <- function(stimuli, inhibitors, signals = NULL, maxStim = 2, maxInhibit = 1) {
##   if (is.null(signals)) { signals <- c(stimuli, inhibitors) }
##   require(CellNOptR)
##   if (maxStim == 0 & maxInhibit == 0) {
##     cnolist <- new("CNOlist",
##                    cues = matrix(0, 1, 1), inhibitors = matrix(0, 1, 1),
##                    stimuli = matrix(0, 1, 1),
##                    signals = list(matrix(0, 1, 1), matrix(0, 1, 1)), timepoints = as.character(c(0,1)))
##   } else {
##     cues <- list()
##     for (i in 1:(length(stimuli)+length(inhibitors))) {
##       cues[[i]] <- c(0,1)
##     }
##     ## expand.grid.jc <- function(x) {
##     ##   temp <- numeric()
##     ##   for (i in 1:(length(x)-1)) {
##     ##     for (j in (i+1):length(x)) {
##     ##       temp <- cbind(temp, cbind(rep.int(x[[i]], length(x[[j]])),
##     ##                                 rep.int(x[[j]], rep.int(length(x[[i]]), length(x[[j]])))))
##     ##     }
##     ##   }
##     ##   return(temp)
##     ## }
##     experiments <- expand.grid(cues)
##     inhibitorSum <- apply(as.matrix(experiments[, (length(stimuli)+1):(length(stimuli) + length(inhibitors))]), 1, sum)
##     experiments <- experiments[which(inhibitorSum <= maxInhibit), ]
##     stimuliSum <- apply(as.matrix(experiments[, 1:length(stimuli)]), 1, sum)
##     experiments <- experiments[which(stimuliSum <= maxStim), ]
##     colnames(experiments)[1:length(stimuli)] <- stimuli
##     colnames(experiments)[(length(stimuli)+1):(length(stimuli)+length(inhibitors))] <- inhibitors
##     experiments <- as.matrix(experiments)
##     rownames(experiments) <- 1:nrow(experiments)
##     for (i in 1:nrow(experiments)) {
##       rownames(experiments)[i] <- paste(colnames(experiments)[which(experiments[i, ] == 1)], collapse = "_")
##     }
##     stimuliTemp <- as.matrix(experiments[, 1:length(stimuli)])
##     colnames(stimuliTemp) <- stimuli
##     experiments2 <- matrix(NA, nrow = nrow(experiments), ncol = length(signals))
##     rownames(experiments2) <- rownames(experiments)
##     colnames(experiments2) <- signals
##     cnolist <- new("CNOlist",
##                    cues = experiments, inhibitors = as.matrix(experiments[, (length(stimuli)+1):(length(stimuli)+length(inhibitors))]),
##                    stimuli = as.matrix(experiments[, 1:length(stimuli)]),
##                    signals = list(experiments2*0, experiments2), timepoints = as.character(c(0,1)))
##   }
##   return(cnolist)
## }
