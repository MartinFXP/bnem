cnoFromData <- function(x, # preprocessed data (no batches allowed!!!)
                         stimuli = c("Tnf", "Bcr"), # stimulated nodes
                         inhibitors = c("Myc", "Ras"), # inhibited genes
                         signals = "none", # "measured genes": if set to "none", inhibitors above are set as inhibited and measured; if set to node names inhibitors are only set as inhibited and signals are set as measured (where egenes can be connected to)
                         times = "none" # c("0h","1h") timepoints
                         ) {
  require(CellNOptR)
  x <- as.matrix(x)
  #x <- x[, order(colnames(x))]
  # make CNOlist:
  CNOlist <- list()
  CNOlist$namesCues <- c(stimuli, inhibitors)
  CNOlist$namesStimuli <- stimuli
  CNOlist$namesInhibitors <- inhibitors
  if (signals[1] == "none") {
    CNOlist$namesSignals <- inhibitors
  } else {
    CNOlist$namesSignals <- signals
  }
  if (times == "none") {
    time <- c(0,1)
    CNOlist$timeSignals <- time
  } else {
    time <- 0:(length(times)-1)
    CNOlist$timeSignals <- time
  }
  for (i in inhibitors) {
    temp <- numeric(ncol(x))
    temp[grep(i, colnames(x))] <- 1
    CNOlist$valueInhibitors <- cbind(CNOlist$valueInhibitors, temp)
  }
  for (i in stimuli) {
    temp <- numeric(ncol(x))
    temp[grep(i, colnames(x))] <- 1
    CNOlist$valueStimuli <- cbind(CNOlist$valueStimuli, temp)
  }
  CNOlist$valueCues <- cbind((CNOlist$valueStimuli), CNOlist$valueInhibitors)
  colnames(CNOlist$valueStimuli) <- stimuli
  colnames(CNOlist$valueInhibitors) <- inhibitors
  colnames(CNOlist$valueCues) <- c(stimuli, inhibitors)
  for (i in 1:length(time)) {
    if (i == 1) {
      ##CNOlist$valueSignals[[as.character(time[i])]] <- matrix(0, nrow = nrow(CNOlist$valueCues), ncol = length(CNOlist$namesSignals))
      CNOlist$valueSignals[[as.character(time[i])]] <- matrix(NA, nrow = nrow(CNOlist$valueCues), ncol = length(CNOlist$namesSignals))
      colnames(CNOlist$valueSignals[[as.character(time[i])]]) <- CNOlist$namesSignals
    } else {
      ##CNOlist$valueSignals[[as.character(time[i])]] <- matrix(sample(c(0,1), nrow(CNOlist$valueCues)*length(CNOlist$namesSignals), replace = TRUE), nrow = nrow(CNOlist$valueCues), ncol = length(CNOlist$namesSignals)) # just random values that are not used later on
      CNOlist$valueSignals[[as.character(time[i])]] <- matrix(NA, nrow = nrow(CNOlist$valueCues), ncol = length(CNOlist$namesSignals)) # just random values that are not used later on
      colnames(CNOlist$valueSignals[[as.character(time[i])]]) <- CNOlist$namesSignals
      rownames(CNOlist$valueSignals[[as.character(time[i])]]) <- colnames(x)
    }
  }
  CNOlist <- new("CNOlist", cues = CNOlist$valueCues, inhibitors = CNOlist$valueInhibitors, stimuli = CNOlist$valueStimuli, signals = CNOlist$valueSignals, timepoints = CNOlist$timeSignals)
  CNOlist <- checkCNOlist(CNOlist)
  return(CNOlist)
}
