# c - s and sk - s were extremly robust against noise (i think)

computeFc <- function (CNOlist, y, test = 1) {
  CompMat <- numeric()
  CompMatNames <- character()
  cnolistStimuli <- apply(CNOlist@stimuli, 1, sum)
  cnolistInhibitors <- apply(CNOlist@inhibitors, 1, sum)
  cnolistCues <- apply(CNOlist@cues, 1, sum)
  maxStim <- max(cnolistStimuli)
  maxKd <- max(cnolistInhibitors)
  grepCtrl <- which(cnolistCues == 0)[1]
  grepStims <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors == 0))
  grepKds <- intersect(which(cnolistStimuli == 0), which(cnolistInhibitors != 0))
  grepStimsKds <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors != 0))
  if (test == 1) {
    ## get ctrl_vs_kd:
    inhibitorsNames <- NULL
    for (i in grepKds) {
      inhibitorsNames <- c(inhibitorsNames, paste(names(which(CNOlist@inhibitors[i, ] >= 1)), collapse = "_"))
    }
    if (length(grepKds) > 0) {
      CompMat <- cbind(CompMat, y[, grepKds] - y[, grepCtrl])
      CompMatNames <- c(CompMatNames, paste("Ctrl_vs_", inhibitorsNames, sep = ""))
    }
    ## get ctrl_vs_stim:
    stimuliNames <- NULL
    for (i in grepStims) {
      stimuliNames <- c(stimuliNames, paste(names(which(CNOlist@stimuli[i, ] >= 1)), collapse = "_"))
    }
    if (length(grepStims) > 0) {
      CompMat <- cbind(CompMat, y[, grepStims] - y[, grepCtrl])
      CompMatNames <- c(CompMatNames, paste("Ctrl_vs_", stimuliNames, sep = ""))
    }
    ## get stim_vs_stim_kd:
    combiNames <- NULL
    for (i in grepStimsKds) {
      combiNames <- c(combiNames, paste(names(which(cbind(CNOlist@stimuli, CNOlist@inhibitors)[i, ] >= 1)), collapse = "_"))
    }
    if (length(grepStimsKds) > 0 & length(grepStims) > 0) {
      CompMat <- cbind(CompMat, y[, rep(grepStimsKds, length(grepStims))] - y[, sort(rep(grepStims, length(grepStimsKds)))])
      orderStims <- order(rep(grepStims, length(grepStimsKds)))
      CompMatNames <- c(CompMatNames, paste(rep(stimuliNames, length(combiNames))[orderStims], "_vs_", rep(combiNames, length(stimuliNames)), sep = ""))
    }
    ## get kd_vs_stim_kd:
    combiNames <- NULL
    for (i in grepStimsKds) {
      combiNames <- c(combiNames, paste(names(which(cbind(CNOlist@inhibitors, CNOlist@stimuli)[i, ] >= 1)), collapse = "_"))
    }
    if (length(grepStimsKds) > 0 & length(grepKds) > 0) {
      CompMat <- cbind(CompMat, y[, rep(grepStimsKds, length(grepKds))] - y[, sort(rep(grepKds, length(grepStimsKds)))])
      orderKds <- order(rep(grepKds, length(grepStimsKds)))
      CompMatNames <- c(CompMatNames, paste(rep(inhibitorsNames, length(combiNames))[orderKds], "_vs_", rep(combiNames, length(inhibitorsNames)), sep = ""))
    }
    colnames(CompMat) <- CompMatNames
    if (sum(duplicated(colnames(CompMat)) == TRUE)) {
      CompMat <- CompMat[, -which(duplicated(colnames(CompMat)) == TRUE)]
    }
  } else {
    CompMat <- numeric()
    CompMatNames <- character()
    cnolistStimuli <- apply(CNOlist@stimuli, 1, sum)
    cnolistInhibitors <- apply(CNOlist@inhibitors, 1, sum)
    cnolistCues <- apply(CNOlist@cues, 1, sum)
    maxStim <- max(cnolistStimuli)
    maxKd <- max(cnolistInhibitors)
    grepCtrl <- which(cnolistCues == 0)[1]
    grepStims <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors == 0))
    grepKds <- intersect(which(cnolistStimuli == 0), which(cnolistInhibitors != 0))
    grepStimsKds <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors != 0))
    ## get ctrl_vs_stim:
    for (i in grepStims) {
      stimNames <- paste(c("Ctrl", "vs", sort(names(which(CNOlist@stimuli[i, ] >= 1)))), collapse = "_")
      if (stimNames %in% CompMatNames) {
        next()
      } else {
        CompMatNames <- c(CompMatNames, stimNames)
        CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
      }
    }
    ## get ctrl_vs_kd:
    for (i in grepKds) {
      kdNames <- paste(c("Ctrl", "vs", sort(colnames(CNOlist@inhibitors)[which(CNOlist@inhibitors[i, ] >= 1)])), collapse = "_")
      if (kdNames %in% CompMatNames) {
        next()
      } else {
        CompMatNames <- c(CompMatNames, kdNames)
        CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
      }
    }
    ## get kd_vs_kd_stim:
    for (i in grepKds) {
      kdNames <- paste(sort(colnames(CNOlist@inhibitors)[which(CNOlist@inhibitors[i, ] >= 1)]), collapse = "_")
      for (j in grepStimsKds) {
        if (paste(sort(colnames(CNOlist@inhibitors)[which(CNOlist@inhibitors[j, ] >= 1)]), collapse = "_") %in% kdNames) {
          stimNames <- paste(c(kdNames, "vs", kdNames, sort(names(which(CNOlist@stimuli[j, ] >= 1)))), collapse = "_")
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
    ##if (type == "model") {
    ## get stim_vs_stim_kd:
    for (i in grepStims) {
      stimNames <- paste(sort(names(which(CNOlist@stimuli[i, ] >= 1))), collapse = "_")
      for (j in grepStimsKds) {
        if (paste(sort(names(which(CNOlist@stimuli[j, ] >= 1))), collapse = "_") %in% stimNames) {
          kdNames <- paste(c(stimNames, "vs", stimNames, sort(colnames(CNOlist@inhibitors)[which(CNOlist@inhibitors[j, ] >= 1)])), collapse = "_")
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
    ##}
    ## if (type == "data") {
    ##   # get (stim - ctrl) - (stimkd - kd): (aka silencing effect)
    ##   for (i in grepStimsKds) {
    ##     silNames <- paste(sort(names(which(CNOlist@cues[i, ] >= 1))), collapse = "_")
    ##     stimNames <- paste(sort(names(which(CNOlist@stimuli[i, ] >= 1))), collapse = "_")
    ##     kdNames <- paste(sort(names(which(CNOlist@inhibitors[i, ] >= 1))), collapse = "_")
    ##     for (j in grepKds) {
    ##       kdNames2 <- paste(sort(names(which(CNOlist@inhibitors[j, ] >= 1))), collapse = "_")
    ##       if (kdNames %in% kdNames2) {
    ##         for (k in grepStims) {
    ##           stimNames2 <- paste(sort(names(which(CNOlist@stimuli[k, ] >= 1))), collapse = "_")
    ##           silNames2 <- c(paste(stimNames2, "vs", stimNames2, kdNames2, sep = "_"))
    ##           if (silNames2 %in% CompMatNames) {
    ##             next()
    ##           }
    ##           if (stimNames %in% stimNames2) {
    ##             silNames2 <- c(paste(stimNames2, "vs", stimNames2, kdNames2, sep = "_"))
    ##             if (silNames2 %in% CompMatNames) {
    ##               next()
    ##             }
    ##             CompMatNames <- c(CompMatNames, silNames2)
    ##             CompMat <- cbind(CompMat, ((y[, j] - y[, grepCtrl]) - (y[, i] - y[, k])))
    ##           }
    ##         }
    ##       }
    ##     }
    ##   }
    ## }
    colnames(CompMat) <- CompMatNames
    CompMat <- CompMat[, sort(colnames(CompMat))]
  }
  return(CompMat)
}
