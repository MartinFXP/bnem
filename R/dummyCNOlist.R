#' creates a general CNOlist object from meta information
#' @param stimuli character vector of stimulated genes
#' @param character vector of inhibited genes
#' @param maxStim maximal number of stimulated genes for a single experiment
#' @param maxInhibit maximal number of inhibited genes for a single experiment
#' @param signals character vector of genes which can directly regulate effect reporters
#' @author Martin Pirkl
#' @return general CNOlist object
#' @export
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signals = c("A", "B","C","D"))
dummyCNOlist <-
function(stimuli = NULL, inhibitors = NULL, maxStim = 0, maxInhibit = 0, signals = NULL) {
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
    if (length(stimuli) == 1) {
      stimDesign <- t(stimDesign)
    }
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
        ## index <- as.vector(stimn:(nrow(combs)+stimn) + (combs - 1)*nrow(stimDesign))
        ## stimDesign[index] <- 1
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
    if (length(inhibitors) == 1) {
      inhibDesign <- t(inhibDesign)
    }
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
        ## index <- as.vector(inhibn:(nrow(combs)+inhibn) + (combs - 1)*nrow(inhibDesign))
        ## inhibDesign[index] <- 1
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
      if (length(inhibitors) == 1) {
        inhibDesign <- t(inhibDesign)
      }
      design <- cbind(stimDesign, inhibDesign)
    }
    if (maxStim == 0 & maxInhibit > 0) {
      stimDesign <- matrix(0, nrow(inhibDesign), stimn)
      if (length(stimuli) == 1) {
        stimDesign <- t(stimDesign)
      }
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
  stimDesign <- as.matrix(design[, 1:ncol(stimDesign)])
  inhibDesign <- as.matrix(design[, (ncol(stimDesign)+1):ncol(design)])
  rownames(design) <- rownames(inhibDesign) <- rownames(stimDesign) <- rownames(signalData) <- c("Ctrl", 2:nrow(design))
  ## for (i in 2:nrow(design)) {
  ##   tmp <- paste(colnames(design)[which(design[i, ] == 1)], collapse = "_")
  ##   rownames(design)[i] <- rownames(inhibDesign)[i] <- rownames(stimDesign)[i] <- rownames(signalData)[i] <- tmp
  ## }
  getRowname <- function(i, M) {
    r <- paste(colnames(M)[which(M[i, ] == 1)], collapse = "_")
    return(r)
  }
  rownames(design)[2:nrow(design)] <- rownames(inhibDesign)[2:nrow(design)] <- rownames(stimDesign)[2:nrow(design)] <- rownames(signalData)[2:nrow(design)] <- unlist(lapply(as.list(2:nrow(design)), getRowname, design))
  if (ncol(stimDesign) == 1) {
    colnames(stimDesign) <- stimuli
  }
  cnolist <- new("CNOlist",
                 cues = design, inhibitors = inhibDesign,
                 stimuli = stimDesign,
                 signals = list(signalData, signalData), timepoints = as.character(c(0,1)))
  cnolist <- checkCNOlist(cnolist)
  return(cnolist)
}
