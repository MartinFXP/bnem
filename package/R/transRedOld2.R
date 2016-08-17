transRedOld2 <-
function(g, max.iter = NULL, verbose = FALSE) { # general transitive reduction:
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    max.iter <- length(v) - 2
  }
  a <- dnf2adj(g)
  g2 <- g#transClose(g, max.iter)
  for (iter in 1:max.iter) {
    for (i in 1:ncol(a)) {
      for (j in 1:ncol(a)) {
        check <- grep(paste(colnames(a)[i], ".*=", colnames(a)[j], sep = ""), g2)
        if (i %in% j | length(check) == 0) { next() }
        for (k in 1:ncol(a)) {
          if (j %in% k | i %in% k) { next() }
          if (length(grep(paste(colnames(a)[i], ".*=", colnames(a)[k], sep = ""), g2)) > 0 & length(grep(paste(colnames(a)[k], ".*=", colnames(a)[j], sep = ""), g2)) > 0) {
            g2 <- g2[-check]
          }
        }
      }
    }
  }
  g3 <- transClose(g2, max.iter)
  if (sum(g %in% g3) > 0) {
    g4 <- g[-which(g %in% g3)]
  }
  g5 <- unique(c(g2, g4))
  return(g5)
}
