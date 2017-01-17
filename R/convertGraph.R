#' converts a disjunctive normal form into a conjunctive normal form and vice versa
#' @param g graph in normal form
#' @author Martin Pirkl
#' @return converted graph normal form
#' @export
#' @examples
#' library(bnem)
#' g <- "A+B=C"
#' g2 <- convertGraph(g)
convertGraph <-
  function(g) { ## input graph as disjunctive normal form like that: c("A+B=D", "C=D", "G+F=U", ...); output is the dual element also in disjunctive normal form;
  g <- sort(g)
  targets <- gsub(".*=", "", g)
  g.new <- NULL
  for (i in unique(targets)) {
    dnf <- list()
    count <- 1
    for (j in g[grep(paste("=", i, sep = ""), g)]) {
      dnf[[count]] <- sort(unique(unlist(strsplit(gsub("=.*", "", j), "\\+"))))
      count <- count + 1
    }
    cnf <- expand.grid(dnf)
    dnf <- NULL
    for (j in 1:dim(cnf)[1]) {
      dnf <- c(dnf, paste(sort(unique(unlist(cnf[j, ]))), collapse = "+"))
    }
    dnf <- paste(sort(dnf), "=", i, sep = "")
    g.new <- c(g.new, dnf)
  }
  vertices <- sort(unique(unlist(strsplit(unlist(strsplit(g.new, "=")), "\\+"))))
  for (i in vertices) {
    if (length(grep(paste(i, ".*", i, ".*=", sep = ""), g.new)) > 0) {
      g.new <- g.new[-grep(paste(i, ".*", i, ".*=", sep = ""), g.new)]
    }
  }
  return(g.new)
}
