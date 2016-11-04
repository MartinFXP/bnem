dataRed <- function(NEMlist, stimuli, stimmult, receptors = NULL, inhibitors, cutoffs, direction = NULL) {

  cutoff <- cutoffs[1]
  cutoff2 <- cutoffs[2]
  cutoff3 <- cutoffs[3]
  cutoff4 <- cutoffs[4]
  
  NEMlist1 <- NEMlist

  if (cutoff > 0) {
    targets <- numeric()
    for (i in 1:length(stimuli)) {
      targets <- c(targets, which(abs(NEMlist1$fc[, paste("Ctrl_vs_", stimuli[i], sep = "")]) > cutoff*stimmult[i]))
    }
    targets1 <- targets
  } else {
    targets1 <- 1:nrow(NEMlist1$fc)
  }
  
  if (abs(cutoff2) > 0 & !is.null(receptors)) {
    targets <- numeric()
    for (j in stimuli) {
      if (is.null(receptors[[j]])) { next() }
      for (i in receptors[[j]]) {
        if (!is.null(direction)) {
          if (direction < 0) {
            targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff2 & NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]*NEMlist1$fc[, paste("Ctrl_vs_", j, sep = "")] < 0))
          }
          if (direction > 0) {
            targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff2 & NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]*NEMlist1$fc[, paste("Ctrl_vs_", j, sep = "")] > 0))
          }
        } else {
          targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff2))
        }
      }
    }
    targets1 <- intersect(targets, targets1)
  } else {
    if (cutoff == 0) {
      targets1 <- 1:nrow(NEMlist1$fc)
    }
  }

  if (abs(cutoff3) > 0) { 
    targets <- numeric()
    for (i in inhibitors) {
      for (j in stimuli) {
        if (length(grep(paste(j, "_vs_", j, "_", i, sep = ""), colnames(NEMlist1$fc))) > 0) {
          if (!is.null(direction)) {
            if (direction < 0) {
              targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff3 & NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]*NEMlist1$fc[, paste("Ctrl_vs_", j, sep = "")] < 0))
            }
            if (direction > 0) {
              targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff3 & NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]*NEMlist1$fc[, paste("Ctrl_vs_", j, sep = "")] > 0))
            }
          } else {
            targets <- c(targets, which(abs(NEMlist1$fc[, paste(j, "_vs_", j, "_", i, sep = "")]) > cutoff3))
          }
        }
      }
    }
    targets1 <- intersect(targets, targets1)
  } else {
    if (cutoff == 0 & cutoff2 == 0) {
      targets1 <- 1:nrow(NEMlist1$fc)
    }
  }
  
  if (abs(cutoff4) > 0) { 
    targets <- numeric()
    for (i in inhibitors) {
      targets <- c(targets, which(abs(NEMlist1$fc[, paste("Ctrl_vs_", i, sep = "")]) > cutoff4))
    }
    targets1 <- intersect(targets, targets1)
  } else {
    if (cutoff == 0 & cutoff2 == 0 & cutoff3 == 0) {
      targets1 <- 1:nrow(NEMlist1$fc)
    }
  }

  NEMlist1$fc <- NEMlist1$fc[unique(targets1), ]
  NEMlist1$exprs <- NEMlist1$exprs[unique(targets1), ]

  if (cutoff == cutoff2 & cutoff2 == cutoff3 & cutoff3 == cutoff4 & cutoff == 0) {
    return(NEMlist)
  } else {
    return(NEMlist1)
  }
}
