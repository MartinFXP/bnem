#' @noRd
graph2adj <-
function(gR) {
    adj.matrix <- matrix(0,
                       length(nodes(gR)),
                       length(nodes(gR))
                       )
    rownames(adj.matrix) <- nodes(gR)
    colnames(adj.matrix) <- nodes(gR)
    for (i in 1:length(nodes(gR))) {
    adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
    }

    return(adj.matrix)
  }
