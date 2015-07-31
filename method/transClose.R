transClose <- function(g, max.iter = NULL) { # general transitive closure:
  if (is.null(max.iter)) {
    max.iter <- length(g)
  }
  for (iter in 1:max.iter) {
    g2 <- NULL
    for (i in g) {
      wtw <- unlist(strsplit(i, "="))
      wtw1 <- unlist(strsplit(wtw[1], "\\+"))
      g.tmp <- g[grep(paste(wtw[2], "=.*|", wtw[2], "\\+*.=", sep = ""), g)]
      for (j in g.tmp) {
        if (length(grep(paste("!", wtw[2], sep = ""), j)) > 0) {
          g2 <- c(g2, gsub(wtw[2], gsub(" ", "", paste("!", wtw1, collapse = "\\+")), gsub(paste("!", wtw[2], sep = ""), wtw[2], j)))
        } else {
          g2 <- c(g2, gsub(wtw[2], wtw[1], j))
        }
      }
    }
    g <- unique(c(g, g2))
  }
  return(g)
}
