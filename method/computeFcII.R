computeFcII <- function (y, stimuli = NULL, inhibitors = NULL, batches, runs, extra = 0) {
  design <- makeDesign(y, stimuli, inhibitors, c(batches, runs))
  CompMat <- numeric()
  CompMatNames <- character()
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
      stimuliSum <- design[, grep(paste(stimuli, collapse = "|"), colnames(design))]
    } else {
      stimuliSum <- apply(design[, grep(paste(stimuli, collapse = "|"), colnames(design))], 1, sum)
    }
  }
  if (is.null(inhibitors) == TRUE) {
    inhibitorsSum <- numeric(nrow(design))
  } else {
    if (length(inhibitors) == 1) {
      inhibitorsSum <- design[, grep(paste(inhibitors, collapse = "|"), colnames(design))]
    } else {
      inhibitorsSum <- apply(design[, grep(paste(inhibitors, collapse = "|"), colnames(design))], 1, sum)
    }
  }
  cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors), collapse = "|"), colnames(design))], 1, sum)
  maxStim <- max(stimuliSum)
  maxKd <- max(inhibitorsSum)
  for (run in runs) {
    for (batch in batches) {
      targetRows <- intersect(which(design[, run] == 1), which(design[, batch] == 1))
      if (length(targetRows) == 0) { next() }
      grepCtrl <- intersect(which(cuesSum == 0), targetRows)[1]
      grepStims <- intersect(intersect(which(stimuliSum != 0), which(inhibitorsSum == 0)), targetRows)
      grepKds <- intersect(intersect(which(stimuliSum == 0), which(inhibitorsSum != 0)), targetRows)
      grepStimsKds <- intersect(intersect(which(stimuliSum != 0), which(inhibitorsSum != 0)), targetRows)
      grepStimsStims <- intersect(intersect(which(stimuliSum >= 2), which(inhibitorsSum == 0)), targetRows)
      # get ctrl_vs_stim:
      for (i in grepStims) {
        if (is.na(grepCtrl)) { next() }
        stimNames <- paste(c("Ctrl", "vs", sort(names(which(stimuliDesign[i, ] >= 1))), batch, run), collapse = "_")
        if (stimNames %in% CompMatNames) {
          next()
        } else {
          CompMatNames <- c(CompMatNames, stimNames)
          CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
        }
      }
      # get ctrl_vs_kd:
      for (i in grepKds) {
        if (is.na(grepCtrl)) { next() }
        kdNames <- paste(c("Ctrl", "vs", sort(names(which(inhibitorsDesign[i, ] >= 1))), batch, run), collapse = "_")
        if (kdNames %in% CompMatNames) {
          next()
        } else {
          CompMatNames <- c(CompMatNames, kdNames)
          CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
        }
      }
      # get kd_vs_kd_stim:
      for (i in grepKds) {
        kdNames <- paste(sort(names(which(inhibitorsDesign[i, ] >= 1))), collapse = "_")
        for (j in grepStimsKds) {
          if (paste(sort(names(which(inhibitorsDesign[j, ] >= 1))), collapse = "_") %in% kdNames) {
            stimNames <- paste(c(kdNames, "vs", kdNames, sort(names(which(stimuliDesign[j, ] >= 1))), batch, run), collapse = "_")
            if (stimNames %in% CompMatNames) {
              next()
            } else {
              CompMatNames <- c(CompMatNames, stimNames)
              CompMat <- cbind(CompMat, (y[, j] - y[, i]))
            }
          } else {
            next()
          }
        }
      }
      # get stim_vs_stim_kd:
      for (i in grepStims) {
        if (is.na(grepCtrl)) { next() }
        stimNames <- paste(sort(names(which(stimuliDesign[i, ] >= 1))), collapse = "_")
        for (j in grepStimsKds) {
          if (paste(sort(names(which(stimuliDesign[j, ] >= 1))), collapse = "_") %in% stimNames) {
            kdNames <- paste(c(stimNames, "vs", stimNames, sort(names(which(inhibitorsDesign[j, ] >= 1))), batch, run), collapse = "_")
            if (kdNames %in% CompMatNames) {
              next()
            } else {
              CompMatNames <- c(CompMatNames, kdNames)
              CompMat <- cbind(CompMat, (y[, j] - y[, i]))
            }
          } else {
            next()
          }
        }
      }
      if (extra == 1) {
        ## get stim_vs_stim_stim:
        for (i in grepStims) {
          if (is.na(grepCtrl)) { next() }
          stimNames <- paste(sort(names(which(stimuliDesign[i, ] >= 1))), collapse = "_")
          for (j in grepStimsStims) {
            if (length(unlist(strsplit(stimNames, "_")) %in% names(which(stimuliDesign[j, ] >= 1))) == length(unlist(strsplit(stimNames, "_"))) & length(unlist(strsplit(stimNames, "_")) %in% names(which(stimuliDesign[j, ] >= 1))) < length(names(which(stimuliDesign[j, ] >= 1)))) {
              stim2Names <- paste(c(stimNames, "vs", sort(names(which(stimuliDesign[j, ] >= 1))), batch, run), collapse = "_")
              if (stim2Names %in% CompMatNames) {
                next()
              } else {
                CompMatNames <- c(CompMatNames, stim2Names)
                CompMat <- cbind(CompMat, (y[, j] - y[, i]))
              }
            } else {
              next()
            }
          }
        }
        ## get ctrl_vs_stim_kd:
        for (i in grepStimsKds) {
          if (is.na(grepCtrl)) { next() }
          stimNames <- paste(sort(names(which(stimuliDesign[i, ] >= 1))), collapse = "_")
          stimkdNames <- paste(c("Ctrl_vs", stimNames, sort(names(which(inhibitorsDesign[i, ] >= 1))), batch, run), collapse = "_")
          if (stimkdNames %in% CompMatNames) {
            next()
          } else {
            CompMatNames <- c(CompMatNames, stimkdNames)
            CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
          }
        }
      }
      # get s - sk - (k - ctrl): not trivial
      ## for (i in stimuli) {
      ##   for (j in inhibitors) {
      ##     name <- 
      ##       CompMatNames <- c(CompMatNames, )
      ##   }
      ## }
    }
  }
  colnames(CompMat) <- CompMatNames
  CompMatCont <- CompMat[, sort(colnames(CompMat))]
  CompMat <- CompMatCont
  return(CompMat)
}
