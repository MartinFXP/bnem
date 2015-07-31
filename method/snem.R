snem <- function(x, stimuli, inhibitors, pkn = NULL, pval = 0.05, cutoff = 0.6, prob = NULL,
                  test = "binom" # more in the future?
                  ) {

  adj.mat <- matrix(0, length(c(stimuli,inhibitors)), length(c(stimuli,inhibitors)))
  diag(adj.mat) <- 0
  colnames(adj.mat) <- c(stimuli,inhibitors)
  rownames(adj.mat) <- c(stimuli,inhibitors)
  AND.mat <- list()

  if (is.null(prob)) { prob <- sum(abs(x) > cutoff)/length(x) }

  for (i in stimuli) {
    
    check1 <- i
    
    AND.mat[[paste(check1, collapse = "_")]] <- adj.mat*0

    if (length(grep(paste("Ctrl_vs_", paste(check1, collapse = "_"), sep = ""), colnames(x))) == 0) { next() }
    
    get <- grep(paste("Ctrl_vs_", paste(check1, collapse = "_"), sep = ""), colnames(x))[1]
    regulated <- which(abs(x[, get]) > cutoff)

    reg.len <- length(regulated)

    for (j in inhibitors) {
      
      check2 <- j

      if (length(grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", paste(check2, collapse = "_"), sep = ""), colnames(x))) == 0) { next() }
      
      get.test <- grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", paste(check2, collapse = "_"), sep = ""), colnames(x))[1]
      regulated.test <- which(abs(x[, get.test]) > cutoff)

      reg.test.len <- length(regulated.test)

      reg.test.int <- intersect(regulated, regulated.test)

      if (reg.test.len == 0) { next() }
      
      if (length(check2) == 1) {

        if (test %in% "binom") {
          test1 <- binom.test(sum(x[reg.test.int, get]*x[reg.test.int, get.test] < 0), reg.test.len, prob, "greater")
        }
        
        if (test1$p.value < pval | sum(x[reg.test.int, get]*x[reg.test.int, get.test] < 0)/reg.test.len == 1) {
          AND.mat[[paste(check1, collapse = "_")]][which(rownames(AND.mat[[paste(check1, collapse = "_")]]) %in% check1), check2] <- 1
        } else {
          next()
        }
      }

      for (k in inhibitors[-which(inhibitors %in% check2)]) {

        get.test2 <- grep(paste(paste(check1, collapse = "_"), "_vs_", paste(check1, collapse = "_"), "_", k, sep = ""), colnames(x))[1]
        regulated.test2 <- which(abs(x[, get.test2]) > cutoff)

        reg.test2.len <- length(regulated.test2)

        reg.test2.int <- intersect(regulated.test2, reg.test.int)
        
        reg.test2.int2 <- intersect(regulated.test2, regulated)

        if (length(reg.test2.int) == 0) { next() }
        
        if (test %in% "binom") {
          test2 <- binom.test(sum(x[reg.test2.int, get.test]*x[reg.test2.int, get.test2] > 0), reg.test2.len, prob, "greater")
          if (sum(x[reg.test2.int2, get]*x[reg.test2.int2, get.test2] < 0) <= sum(x[reg.test.int, get]*x[reg.test.int, get.test] < 0)) {
            test3 <- list()
            test3$p.value <- 0
          } else {
            test3 <- binom.test(sum(x[reg.test2.int, get.test]*x[reg.test2.int, get.test2] > 0), sum(x[reg.test.int, get.test]*x[reg.test.int, get.test2] > 0), p = 0.95, alternative = "greater")
          }
        }
        
        if ((test2$p.value < pval & test3$p.value < pval) | sum(x[reg.test2.int, get.test]*x[reg.test2.int, get.test2] > 0)/reg.test2.len == 1) {
          AND.mat[[paste(check1, collapse = "_")]][check2, k] <- 1
        }
      }
    }
  }

  AND.mat.merge <- adj.mat*0

  for (i in 1:length(AND.mat)) {

    AND.mat.merge <- AND.mat.merge + AND.mat[[i]]

  }

  AND.mat.merge[which(AND.mat.merge > 1)] <- 1

  if (!is.null(pkn)) {
    AND.mat.merge <- AND.mat.merge*pkn
  }
  
  AND.mat.merge <- transitive.closure(AND.mat.merge, mat = T)

  return(AND.mat.merge)

}
