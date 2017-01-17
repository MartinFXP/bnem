#' @noRd
drawLocal <-
function(l, type = "l", ...) {
  local.edges <- c(1, l$edges[[1]])
  local.edges1 <- local.edges
  if (length(unique(local.edges)) < 2) {
    local.edges <- rep(mean(l$scores[[1]]), length(local.edges1))
    min.local.edges <- min(local.edges)
    convert.edges <- seq(min(l$scores[[1]]), max(l$scores[[1]]), length.out = max(local.edges - min.local.edges))
    convert.legend <- sort(unique(local.edges))
  } else {
    min.local.edges <- min(local.edges)
    convert.edges <- seq(min(l$scores[[1]]), max(l$scores[[1]]), length.out = max(local.edges - min.local.edges)+1)
    convert.legend <- sort(unique(local.edges))
    for (i in 1:(max(local.edges - min.local.edges)+1)) {
      local.edges[which(local.edges == (i+min.local.edges)-1)] <- convert.edges[which(convert.legend == (i+min.local.edges)-1)]
    }
  }
  ylim <- c(min(l$scores[[1]]), max(l$scores[[1]]))
  par(mar=c(4, 4, 4, 4) + 0.1)
  plot(l$scores[[1]], type = type, main = "score improvement (black) and number of changed edges (red)", ylab = "score", xlab = "move", ...)
  lines(local.edges, col = 2, type = "b")
  axis(4, at = sort(unique(local.edges)), labels = sort(unique(local.edges1)), col = "red", col.ticks = "red", col.axis = "red")
  mtext("edges changed", side=4, line=3, cex.lab=1,las=0, col="red")
  ##axis(4, at = sort(unique(local.edges)), labels = sort(unique(local.edges1)))
  #abline(h = convert.edges[which(convert.legend == 1)], col = 2)
}
