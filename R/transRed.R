#' calculates transitive reduciton of a hyper-graph in normal form
#' @param g hyper-graph in normal form
#' @param max.iter maximal number of iterations till convergence
#' @param verbose verbose output?
#' @author Martin Pirkl
#' @return transitive reduction of the hyper-graph in normal form
#' @export
#' @examples
#' g <- c("A=B", "A=C", "B=C", "B=D", "!A=D")
#' gred <- transRed(g)
transRed <-
function(g, max.iter = NULL, verbose = FALSE) { # general transitive reduction:
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    max.iter <- length(v) - 2
  }
  a <- dnf2adj(g)
  g2 <- g
  h <- getHierarchy(g2)
  if (length(h) > 2) {
    for (i in 1:(length(h)-2)) {
      for (j in h[[i]]) {
        for (k in (i+2):length(h)) {
          for (l in h[[k]]) {
            if (length(grep(paste(".*", j, ".*=", l, sep = ""), g2)) != 0) {
              if (length(grep(paste(".*=", l, sep = ""), g2)) > 1) {
                g2 <- g2[-grep(paste(".*", j, ".*=", l, sep = ""), g2)]
              }
            }
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
