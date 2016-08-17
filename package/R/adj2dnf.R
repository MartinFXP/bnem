adj2dnf <-
function(A) {

  dnf <- NULL
  
  for (i in 1:ncol(A)) {
    for (j in 1:nrow(A)) {
      if (i %in% j) { next() }
      if (A[i, j] == 1) {
        dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
      }
      if (A[i, j] == -1) {
        dnf <- c(dnf, paste("!", colnames(A)[i], "=", rownames(A)[j], sep = ""))
      }
    }
  }

  dnf <- unique(dnf)
  
  return(dnf)

}
