#' @noRd
disc <-
function(x, beta = 0.5) { # simple discretization
  x.disc <- x
  x.disc[which(abs(x.disc) <= beta)] <- 0
  x.disc[which(x.disc > beta)] <- 1
  x.disc[which(x.disc < -beta)] <- -1
  return(x.disc)
}
