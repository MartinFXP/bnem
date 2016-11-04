sbnem <- function(x, stimuli, inhibitors, pkn = NULL, prob = NULL, cutoff = 0.6, pval = 0.05) {

  adj.mat <- matrix(0, length(c(stimuli,inhibitors)), length(c(stimuli,inhibitors)))
  diag(adj.mat) <- 0
  colnames(adj.mat) <- c(stimuli,inhibitors)
  rownames(adj.mat) <- c(stimuli,inhibitors)
  AND.mat <- list()

  if (is.null(prob)) { prob <- sum(abs(x) > cutoff)/length(x) }

  for (i in stimuli)  {#1:bin2dec(rep(1,length(stimuli)))) {
    #binary1 <- dec2bin(i)
    #binary1 <- c(rep(0, length(stimuli) - length(binary1)), binary1)
    #binary1 <- binary1[length(binary1):1]
    check1 <- i#stimuli[which(binary1 == 1)]
    
    AND.mat[[paste(check1, collapse = "_")]] <- adj.mat*0

    if (length(grep(paste("Ctrl_vs_", paste(check1, collapse = "_"), sep = ""), colnames(x))) == 0) { next() }
    
    stim.col <- grep(paste("Ctrl_vs_", paste(check1, collapse = "_"), sep = ""), colnames(x))[1] # should only be one anyway
    
    stim.targets <- which(abs(x[, stim.col]) > cutoff)

    for (j in inhibitors) {#1:bin2dec(rep(1,length(inhibitors)))) {

                                        #binary2 <- dec2bin(j)
                                        #binary2 <- c(rep(0, length(inhibitors) - length(binary2)), binary2)
                                        #binary2 <- binary2[length(binary2):1]
      check2 <- j#inhibitors[which(binary2 == 1)]

      if (length(grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", paste(check2, collapse = "_"), sep = ""), colnames(x))) == 0) { next() }
      
      gene.col <- grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", paste(check2, collapse = "_"), sep = ""), colnames(x))[1]

      gene.targets <- which(x[which(abs(x[, gene.col]) > cutoff), gene.col]*x[which(abs(x[, gene.col]) > cutoff), stim.col] < 0)

      if (length(gene.targets)/sum(abs(x[, gene.col]) > cutoff) >= prob) {
        AND.mat[[paste(check1, collapse = "_")]][check1, check2] <- 1
      }
      
      for (k in inhibitors[-which(inhibitors %in% check2)]) {
        
        if (length(grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", k, sep = ""), colnames(x))) == 0) { next() }

        gene2.col <- grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", k, sep = ""), colnames(x))[1]

        gene2.targets <- which(x[which(abs(x[, gene2.col]) > cutoff), gene2.col]*x[which(abs(x[, gene2.col]) > cutoff), stim.col] < 0)
        overlap <- intersect(names(gene.targets), names(gene2.targets))

        gene.gene2 <- which(x[which(rownames(x) %in% overlap), gene.col]*x[which(rownames(x) %in% overlap), gene2.col] > 0)

        if (binom.test(length(gene.gene2), length(gene2.targets), p = prob, alternative = "greater")$p.value < 0.05) {

          AND.mat[[paste(check1, collapse = "_")]][check2, k] <- 1

        }
        
      }

    }

  }

  return(AND.mat[[1]])

}
