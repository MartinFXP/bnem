addSignal <-
function(s, CNOlist, stim = NULL, inhibit = NULL) {
  CNOlist2 <- CNOlist
  CNOlist2@signals[[1]] <- cbind(CNOlist2@signals[[1]], s = 0)
  CNOlist2@signals[[2]] <- cbind(CNOlist2@signals[[2]], s = 0)
  colnames(CNOlist2@signals[[1]]) <- c(colnames(CNOlist@signals[[1]]), s)
  colnames(CNOlist2@signals[[2]]) <- c(colnames(CNOlist@signals[[2]]), s)
  if (!is.null(inhibit)) {
    CNOlist2@cues <- cbind(CNOlist2@cues, 0)
    CNOlist2@cues <- rbind(CNOlist2@cues, matrix(0, nrow = nrow(inhibit), ncol = ncol(CNOlist2@cues)))
    CNOlist2@cues[(nrow(CNOlist2@cues) - nrow(inhibit) + 1):nrow(CNOlist2@cues), which(colnames(CNOlist2@cues) %in% colnames(inhibit))] <- inhibit
    CNOlist2@cues[(nrow(CNOlist2@cues) - nrow(inhibit) + 1):nrow(CNOlist2@cues), ncol(CNOlist2@cues)] <- 1
    colnames(CNOlist2@cues)[ncol(CNOlist2@cues)] <- s
    CNOlist2@stimuli <- CNOlist2@cues[, which(colnames(CNOlist2@cues) %in% colnames(CNOlist2@stimuli))]
    CNOlist2@inhibitors <- CNOlist2@cues[, which(colnames(CNOlist2@cues) %in% c(colnames(CNOlist2@inhibitors), s))]
    CNOlist2@signals[[1]] <- rbind(CNOlist2@signals[[1]], matrix(0, nrow = nrow(inhibit), ncol = ncol(CNOlist2@signals[[1]])))
    CNOlist2@signals[[2]] <- rbind(CNOlist2@signals[[2]], matrix(0, nrow = nrow(inhibit), ncol = ncol(CNOlist2@signals[[2]])))
  }
  CNOlist2 <- checkCNOlist(CNOlist2)
  return(CNOlist2)
}
