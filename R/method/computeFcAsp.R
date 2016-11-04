# c - s and sk - s were extremly robust against noise (i think)

computeFcAsp <- function (CNOlist,
                       y,
                       file = "dataFc.lp") {
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
  # get ctrl_vs_stim:
  for (i in grepStims) {
    stimNames <- paste(c("Ctrl", "vs", sort(names(which(CNOlist@stimuli[i, ] >= 1)))), collapse = "_")
    if (stimNames %in% CompMatNames) {
      next()
    } else {
      CompMatNames <- c(CompMatNames, stimNames)
      CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
      write(paste("datafc(", grepCtrl, ",", i, ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(y)))))), ",", disc((y[, i] - y[, grepCtrl])), ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
      write(paste("comp(", grepCtrl, ",", i, ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
    }
  }
  # get ctrl_vs_kd:
  for (i in grepKds) {
    kdNames <- paste(c("Ctrl", "vs", sort(names(which(CNOlist@inhibitors[i, ] >= 1)))), collapse = "_")
    if (kdNames %in% CompMatNames) {
      next()
    } else {
      CompMatNames <- c(CompMatNames, kdNames)
      CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
      write(paste("datafc(", grepCtrl, ",", i, ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(y)))))), ",", disc((y[, i] - y[, grepCtrl])), ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
      write(paste("comp(", grepCtrl, ",", i, ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
    }
  }
  # get kd_vs_kd_stim:
  for (i in grepKds) {
    kdNames <- paste(sort(names(which(CNOlist@inhibitors[i, ] >= 1))), collapse = "_")
    for (j in grepStimsKds) {
      if (paste(sort(names(which(CNOlist@inhibitors[j, ] >= 1))), collapse = "_") %in% kdNames) {
        stimNames <- paste(c(kdNames, "vs", kdNames, sort(names(which(CNOlist@stimuli[j, ] >= 1)))), collapse = "_")
        if (stimNames %in% CompMatNames) {
          next()
        } else {
          CompMatNames <- c(CompMatNames, stimNames)
          CompMat <- cbind(CompMat, (y[, j] - y[, i]))
          write(paste("datafc(", i, ",", j, ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(y)))))), ",", disc((y[, j] - y[, i])), ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
          write(paste("comp(", i, ",", j, ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
        }
      } else {
        next()
      }
    }
  }
  #if (type == "model") {
    # get stim_vs_stim_kd:
    for (i in grepStims) {
      stimNames <- paste(sort(names(which(CNOlist@stimuli[i, ] >= 1))), collapse = "_")
      for (j in grepStimsKds) {
        if (paste(sort(names(which(CNOlist@stimuli[j, ] >= 1))), collapse = "_") %in% stimNames) {
          kdNames <- paste(c(stimNames, "vs", stimNames, sort(names(which(CNOlist@inhibitors[j, ] >= 1)))), collapse = "_")
          if (kdNames %in% CompMatNames) {
            next()
          } else {
            CompMatNames <- c(CompMatNames, kdNames)
            CompMat <- cbind(CompMat, (y[, j] - y[, i]))
            write(paste("datafc(", i, ",", j, ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(y)))))), ",", disc((y[, j] - y[, i])), ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
            write(paste("comp(", i, ",", j, ").", sep = ""), file = paste("dissertation/asp_implementation/dataFc.lp", sep = ""), append = TRUE)
          }
        } else {
          next()
        }
      }
    }
  #}
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
  return(CompMat)
}
