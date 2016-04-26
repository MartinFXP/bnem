transClose <- function(g, max.iter = NULL, verbose = FALSE) { # does soemthign strange!!!
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    max.iter <- length(v) - 2
  }
  if (verbose) {
    print(paste("maximum iterations: ", max.iter, sep = ""))
  }
  g.out <- unique(gsub(".*=", "", g))
  g.closed <- g
  for (iter in 1:max.iter) {
    g.old <- g.closed
    
    if (verbose) {
      cat('\r', paste("iteration: ", iter, sep = ""))
      flush.console()
    }
    for (i in g.closed) {
      input <- gsub("!", "", unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+")))
      input <- intersect(input, g.out)
      output <- gsub(".*=", "", i)
      if (length(input) == 0) { next() }
      for (j in unique(input)) {
        if (j %in% unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+"))) {
          for (k in gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)])) {
            g.closed <- c(g.closed, gsub(j, k, i))
          }
        } else {
          literals <- list()
          count <- 0
          for (k in unique(gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)]))) {
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
      if (verbose) {
        cat("\n")
        print(paste("successfull convergence", sep = ""))
      }
      break()
    }
  }
  if (verbose) {
    cat("\n")
  }
  g.closed <- unique(g.closed)
  for (j in 1:length(g.closed)) {
    i <- g.closed[j]
    input <- unlist(strsplit(i, "="))
    output <- input[2]
    input <- unlist(strsplit(input[1], "\\+"))
    input <- unique(input)
    g.closed[j] <- paste(paste(input, collapse = "+"), output, sep = "=")
  }
  if (length(grep(paste(paste(v, ".*", v, ".*=", sep = ""), collapse = "|"), g.closed)) > 0) {
    g.closed <- g.closed[-grep(paste(paste(v, ".*", v, ".*=", sep = ""), collapse = "|"), g.closed)]
  }
  return(g.closed)
}
