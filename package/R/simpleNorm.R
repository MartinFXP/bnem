simpleNorm <-
function(x) {
  if (is.matrix(x) == FALSE) {
    genes.min <- min(x)
    x <- (x - genes.min)
    genes.max <- max(x)
    x <- x/genes.max
  } else {
    genes.min <- apply(x, 1, min)
    x <- (x - genes.min)
    genes.max <- apply(x, 1, max)
    x <- x/genes.max
  }
  x[is.na(x)] <- 0
  return(x)
}
