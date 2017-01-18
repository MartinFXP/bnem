#' Convert epiNEM model into general Boolean graph.
#' Only needed for comparing accuracy of inferred network for bnem and epiNEM.
#' @param t full epiNEM model
#' @author Martin Pirkl
#' @seealso CreateTopology
#' @export
#' @examples
#' library(epiNEM)
#' topology <- CreateTopology(3, 1, force = TRUE)
#' topology <- unlist(unique(topology), recursive = FALSE)
#' extTopology <- ExtendTopology(topology$model, 100)
#' b <- EpiNEM2BooleanGraph(extTopology)
#' @return boolean hyper-graph
epiNEM2Bg <- function(t) {
    if (is.matrix(t)) {
        colnames(t) <- paste("S_vs_S_", gsub("\\.", "_", colnames(t)), sep = "")
        return(t)
    } else {
        tmp <- apply(t$origModel, 2, sum)
        stim <- rownames(t$origModel)[which(tmp == min(tmp))]
        graph <- NULL

        for (i in 1:length(t$column)) {
            parents <- sort(rownames(t$origModel)[which(t$origModel[, t$column[i]] == 1)])
            child <- colnames(t$origModel)[t$column[i]]
            if (length(parents) == 2) {
                if (t$logics[i] %in% "OR") {
                    graph <- unique(c(graph, convertGraph(adj2dnf(t$origModel))))
                }
                if (t$logics[i] %in% "AND") {
                    graph <- unique(c(graph, transRed(convertGraph(adj2dnf(t$origModel)))))
                    graph <- c(graph, convertGraph(graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))
                    graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]
                }
                if (t$logics[i] %in% paste(parents[2], " masks the effect of ", parents[1], sep = "")) {
                    graph <- c(graph, unique(convertGraph(adj2dnf(t$origModel))))
                    graph <- c(graph, gsub(parents[2], paste("!", parents[2], sep = ""), gsub("\\+\\+|^\\+", "", gsub(parents[1], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))), gsub("\\+=", "=", gsub("\\+\\+|^\\+", "", gsub(parents[2], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))))
                    graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]
                }
                if (t$logics[i] %in% paste(parents[1], " masks the effect of ", parents[2], sep = "")) {
                    graph <- c(graph, unique(convertGraph(adj2dnf(t$origModel))))
                    graph <- c(graph, gsub(parents[1], paste("!", parents[1], sep = ""), gsub("\\+\\+|^\\+", "", gsub(parents[2], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))), gsub("\\+=", "=", gsub("\\+\\+|^\\+", "", gsub(parents[1], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))))
                    graph <- gsub("\\+=", "=", graph)
                    graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]
                }
                if (t$logics[i] %in% "XOR") {
                    graph <- unique(c(graph, convertGraph(adj2dnf(t$origModel))))
                    edge <-  graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]
                    for (j in parents) {
                        edge <- gsub(j, paste("!", j, sep = ""), edge)
                    }
                    graph <- unique(c(graph, edge))
                }
            }
            if (length(parents) > 2 | length(parents) == 1) {
                graph <- c(graph, paste(sort(parents), child, sep = "="))
            }
        }
        all <- rownames(t$origModel)
        children2 <- unique(gsub(".*=", "", graph))
        if (sum(!(all %in% children2)) > 0) {
            graph <- c(graph, paste("S", all[which(!(all %in% children2))], sep = "="))
        }
        return(unique(graph))
    }
}
