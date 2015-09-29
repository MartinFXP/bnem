transClose <- function(g, max.iter = NULL) {
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    max.iter <- length(v) - 2
  }
  print(paste("maximum iterations: ", max.iter, sep = ""))
  g.out <- unique(gsub(".*=", "", g))
  g.closed <- g
  for (iter in 1:max.iter) {
    g.old <- g.closed
    
    cat('\r', paste("iteration: ", iter, sep = ""))
    flush.console()
        
    for (i in g.closed) {
      input <- gsub("!", "", unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+")))
      input <- intersect(input, g.out)
      output <- gsub(".*=", "", i)
      if (length(input) == 0) { next() }
      for (j in input) {
        if (j %in% unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+"))) {
          tmp <- paste(sort(unique(unlist(strsplit(gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)]), "\\+")))), collapse = "+")
          g.closed <- c(g.closed, gsub(j, tmp, i))
        } else {
          literals <- list()
          count <- 0
          for (k in gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)])) {
            count <- count + 1
            literals[[count]] <- gsub("!!", "", paste("!", unlist(strsplit(k, "\\+")), sep = ""))
          }
          combis <- expand.grid(literals)
          combis <- apply(combis, c(1,2), as.character)
          for (k in 1:nrow(combis)) {
            g.closed <- c(g.closed, gsub(paste("!", j, sep = ""), paste(combis[k, ], collapse = "+"), i))
          }
        }
      }
    }
    if (all(g.closed %in% g.old)) {
      cat("\n")
      print(paste("successfull convergence", sep = ""))
      break()
    }
  }
  cat("\n")
  g.closed <- unique(g.closed)
  return(g.closed)
}
