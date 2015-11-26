randomDnf <- function(vertices = 10, negation = TRUE, max.edge.size = NULL, max.edges = NULL, dag = FALSE) {
  dnf <- NULL
  if (is.numeric(vertices)) {
    vertices <- LETTERS[1:vertices]
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
    dnf <- c(dnf, paste(c(paste(paste(pre, inputs, sep = ""), collapse = "+"), "=", output), collapse = ""))
  }
  if (dag) {
    dnf <- removeCycles(dnf = dnf)
  }
  return(dnf)
}

## dnf <- randomDnf(5, dag = T); par(mfrow=c(1,3)); plotDnf(dnf); plotDnf(transRed(dnf)); plotDnf(transClose(dnf));
