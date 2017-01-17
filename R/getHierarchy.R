#' @noRd
getHierarchy <-
function(graph) {
  adj <- dnf2adj(graph)
  dnf <- adj2dnf(adj)
  g <- plotDnf(dnf, draw = FALSE)
  Ypos <- g@renderInfo@nodes$labelY
  Ynames <- names(g@renderInfo@nodes$labelY)
  ## Ypos <- Ypos[-grep("and", Ynames)]
  ## Ynames <- Ynames[-grep("and", Ynames)]
  hierarchy <- list()
  count <- 0
  for (i in sort(unique(Ypos), decreasing = TRUE)) {
    count <- count + 1
    hierarchy[[count]] <- Ynames[which(Ypos == i)]
  }
  return(hierarchy)
}
