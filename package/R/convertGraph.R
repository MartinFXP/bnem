convertGraph <-
function(g) { ## input graph as disjunctive normal form like that: c("A+B=D", "C=D", "G+F=U", ...); output is the dual element also in disjunctive normal form;
  targets <- gsub(".*=", "", g)
  g.new <- NULL
  for (i in unique(targets)) {
    dnf <- list()
    count <- 1
    for (j in g[grep(paste("=", i, sep = ""), g)]) {
      dnf[[count]] <- unlist(strsplit(gsub("=.*", "", j), "\\+"))
      count <- count + 1
    }
    cnf <- expand.grid(dnf)
    dnf <- NULL
    for (j in 1:dim(cnf)[1]) {
      dnf <- c(dnf, paste(unique(unlist(cnf[j, ])), collapse = "+"))
    }
    dnf <- paste(dnf, "=", i, sep = "")
    g.new <- c(g.new, dnf)
  }
  vertices <- unique(unlist(strsplit(unlist(strsplit(g.new, "=")), "\\+")))
  for (i in vertices) {
    if (length(grep(paste(i, ".*", i, ".*=", sep = ""), g.new)) > 0) {
      g.new <- g.new[-grep(paste(i, ".*", i, ".*=", sep = ""), g.new)]
    }
  }
  return(g.new)
}
