#' creates a random normal form or hyper-graph
#' @param maximal number of vertices
#' @param negation allowed?
#' @param max.edge.size maximal number of inputs per edge
#' @param max.edges maximal number of hyper-edges
#' @param dag is the graph to ba a dag?
#' @author Martin Pirkl
#' @return random hyper-graph in normal form
#' @export
#' @examples
#' library(bnem)
#' g <- randomDnf(10)
randomDnf <-
function(vertices = 10, negation = TRUE, max.edge.size = NULL, max.edges = NULL, dag = FALSE) {
  dnf <- NULL
  if (is.numeric(vertices)) {
    if (vertices < 27) {
      vertices <- LETTERS[1:vertices]
    } else {
      vertices <- paste("S", 1:vertices, "g", sep = "")
    }
  }
  if (is.null(max.edge.size)) {
    max.edge.size <- length(vertices) - 1
  }
  if (is.null(max.edges)) {
    max.edges <- length(vertices) - 1
  }
  for (i in 1:max.edges) {
    edge.size <- sample(1:max.edge.size, 1)
    output <- sample(vertices, 1)
    inputs <- NULL
    for (j in 1:edge.size) {
      inputs <- c(inputs, sample(vertices[-grep(paste(c(output, inputs), collapse = "|"), vertices)], 1))
    }
    if (negation) {
      pre <- sample(c("", "!"), edge.size, replace = T)
    } else {
      pre <- rep("", edge.size)
    }
    inputs <- sort(inputs)
    dnf <- c(dnf, paste(c(paste(paste(pre, inputs, sep = ""), collapse = "+"), "=", output), collapse = ""))
  }
  if (dag) {
    dnf <- removeCycles(dnf = dnf)
  }
  return(unique(dnf))
}
