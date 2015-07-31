plotDnf <- function(dnf = NULL, freq = NULL, stimuli = c(), signals = c(), inhibitors = c(), connected = TRUE,  CNOlist = NULL, cex = NULL, fontsize = NULL, labelsize = NULL, type = 2, lwd = 2, edgelwd = 2, legend = 0, x = 0, y = 0, xjust = 0, yjust = 0, width = 1.5, height = 1, rankdir = "TB", rank = "same", layout = "dot", main = "", sub = "", cex.main = 1.5, cex.sub = 1, col.sub = "grey", fontcolor = NULL, nodestates = NULL, simulate = NULL, andcolor = "transparent", edgecol = NULL, labels = NULL, labelcol = "blue", nodelabel = NULL, nodecol = NULL, bordercol = NULL, nodeshape = NULL, verbose = FALSE, edgestyle = NULL, nodeheight = NULL, nodewidth = NULL, edgewidth = NULL, lty = NULL, hierarchy = NULL, nodefontsize = NULL, edgehead = NULL, edgelabel = NULL, ...) {
  ## see graphvizCapabilities()$layoutTypes for supported layouts
  require(Rgraphviz)

  graph <- dnf

  if (!is.null(hierarchy)) {
    graph2 <- NULL
    for (i in graph) {
      input <- unlist(strsplit(i, "="))
      output <- input[2]
      input <- gsub("!", "", unlist(strsplit(input[1], "\\+")))
      for (j in input) {
        graph2 <- c(graph2, paste(j, output, sep = "="))
      }
      graph2 <- unique(graph2)
    }
    hgraph <- NULL
    for (i in 1:(length(hierarchy)-1)) {
      for (j in hierarchy[[i]]) {
        for (k in hierarchy[[i+1]]) {
          hgraph <- c(hgraph, paste(j, k, sep = "="))
        }
      }
    }
    if (sum(hgraph %in% graph2) > 0) {
      hgraph <- hgraph[-which(hgraph %in% graph2)]
    }
    dnf <- c(hgraph, graph)
    ## update all the parameters e.g. edgestyle...
    if (is.null(edgecol)) {
      edgecol <- c(rep("transparent", length(hgraph)), rep("black", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
    } else {
      if (length(edgecol) == 1) {
        edgecol <- c(rep("transparent", length(hgraph)), rep(edgecol, length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
      } else {
        edgecol <- c(rep("transparent", length(hgraph)), edgecol)
      }
    }
    if (is.null(lty)) {
      lty <- c(rep("solid", length(hgraph)), rep("solid", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
    } else {
      lty <- c(rep("solid", length(hgraph)), lty)
    }
    if (is.null(edgewidth)) {
      edgewidth <- c(rep(1, length(hgraph)), rep(1, length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
    } else {
      edgewidth <- c(rep(1, length(hgraph)), edgewidth)
    }
    if (!is.null(labels)) {
      labels <- c(rep(1, length(hgraph)), labels)
    }
  } else {
    if (is.null(lty) & !is.null(dnf)) {
      lty <- c(rep("solid", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
    }
  }

  dolegend <- FALSE
  if (is.null(dnf)) {
    dnf <- c("A=B")
    dolegend <- TRUE
  }
  
  if (!is.null(simulate)) {
    nodestates <- simulateDnf(graph, stimuli = simulate$stimuli, inhibitors = simulate$inhibitors)
  }
  
  if (is.null(freq)) {
    use.freq = FALSE
  } else {
    use.freq = TRUE
    if (is.null(labels)) {
      labels <- as.character(round(freq, 2)*100)
    }
  }
  if (is.null(labels)) {
    labels <- rep("", length(dnf))
  }
  
  if (is.null(fontsize)) {
    fontsize <- ""
  }
  if (is.null(labelsize)) {
    labelsize <- fontsize
  }

  if (!is.null(CNOlist)) {
    if (length(stimuli) == 0) {
      stimuli <- colnames(CNOlist@stimuli)
    }
    if (length(signals) == 0) {
      signals <- colnames(CNOlist@signals[[1]])
    }
    if(length(inhibitors) == 0) {
      inhibitors <- colnames(CNOlist@inhibitors)
    }
  }

  if (connected) {
    Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))))
    if (length(grep("\\+", dnf)) > 0) {
      Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
    }
    V <- unique(gsub("!", "", Vneg))
    stimuli <- intersect(stimuli, V)
    signals <- intersect(signals, V)
    inhibitors <- intersect(inhibitors, V)
  } else {
    Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))), stimuli, signals, inhibitors))
    if (length(grep("\\+", dnf)) > 0) {
      Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
    }
    V <- unique(gsub("!", "", Vneg))
  }

  V <- sort(V)
  
  Vneg <- c(V, Vneg[grep("!", Vneg)])

  if (!is.null(nodecol)) {
    if (length(nodecol) == 1 & !is.list(nodecol)) {
      nodecol.tmp <- nodecol
      nodecol <- list()
      for (i in V[-grep("and", V)]) {
        nodecol[[i]] <- nodecol.tmp
      }
    }
  }

  if (!is.null(bordercol)) {
    if (length(bordercol) == 1 & !is.list(bordercol)) {
      bordercol.tmp <- bordercol
      bordercol <- list()
      for (i in V) {
        if (length(grep("and", i)) == 0) {
          bordercol[[i]] <- bordercol.tmp
        }
      }
    }
  }

  E <- list()

  for (i in V) {
    E[[i]] <- list()
  }

  Eneg <- list()

  for (i in Vneg) {
    Eneg[[i]] <- list()
  }

  count <- 0
  
  for (i in dnf) {
    print(i)
    tmp <- unlist(strsplit(i, "="))
    if (length(tmp)==1) {
      Eneg[[tmp]][["edges"]] <- c(Eneg[[tmp]][["edges"]], NULL)
      tmp <- gsub("!", "", tmp)
      E[[tmp]][["edges"]] <- c(E[[tmp]][["edges"]], NULL)
    } else {
      tmp2 <- unlist(strsplit(tmp[1], "\\+"))
      if (length(tmp2) > 1) {
        count <- count + 1
        Eneg[[paste("and", count, sep = "")]][["edges"]] <- c(Eneg[[paste("and", count, sep = "")]][["edges"]], which(Vneg %in% tmp[2]))
        E[[paste("and", count, sep = "")]][["edges"]] <- c(E[[paste("and", count, sep = "")]][["edges"]], which(V %in% tmp[2]))
        for (j in tmp2) {
          Eneg[[j]][["edges"]] <- c(Eneg[[j]][["edges"]], which(Vneg %in% paste("and", count, sep = "")))
          j <- gsub("!", "", j)
          E[[j]][["edges"]] <- c(E[[j]][["edges"]], which(V %in% paste("and", count, sep = "")))
        }
      } else {
        print(tmp)
        print(tmp2)
        Eneg[[tmp2]][["edges"]] <- c(Eneg[[tmp2]][["edges"]], which(Vneg %in% tmp[2]))
        tmp2 <- gsub("!", "", tmp2)
        E[[tmp2]][["edges"]] <- c(E[[tmp2]][["edges"]], which(V %in% tmp[2]))
      }
    }
  }

  print(names(Eneg))

  g <- new("graphNEL",nodes=V,edgeL=E,edgemode="directed")

  gneg <- new("graphNEL",nodes=Vneg,edgeL=Eneg,edgemode="directed")

  nodes <- buildNodeList(g)

  edges <- buildEdgeList(g)
  
  nodesneg <- buildNodeList(gneg)

  edgesneg <- buildEdgeList(gneg)

  edgesnew <- list()
  
  for (i in sort(names(edges))) {
    edgesnew <- c(edgesnew, edges[[i]])
  }

  names(edgesnew) <- sort(names(edges))

  edges <- edgesnew

  edgesnegnew <- list()
  
  for (i in names(edgesneg)[order(gsub("!", "", names(edgesneg)))]) {
    edgesnegnew <- c(edgesnegnew, edgesneg[[i]])
  }

  names(edgesnegnew) <- names(edgesneg)[order(gsub("!", "", names(edgesneg)))]

  edgesneg <- edgesnegnew

  print(names(edges))
  print(names(edgesneg))

  if (verbose) {
    print(paste("order of nodes: ", paste(names(nodes), collapse = ", "), sep = ""))
    print(paste("order of edges: ", paste(names(edges), collapse = ", "), sep = ""))
  }

  nodeshape2 <- nodeshape
  nodeshape <- list()
  if (length(nodeshape2) == 1 & !(is.list(nodeshape2))) {
    for (i in 1:length(nodes)) {
      nodeshape[[nodes[[i]]@name]] <- nodeshape2
    }
  } else {
    nodeshape <- nodeshape2
  }
  
  for (i in 1:length(nodes)) {
    nodes[[i]]@attrs$height <- height
    nodes[[i]]@attrs$width <- width
    if (!is.null(nodelabel[[nodes[[i]]@name]])) {
      nodes[[i]]@attrs$name <- nodelabel[[nodes[[i]]@name]]
    }
    if (!is.null(nodeheight[[nodes[[i]]@name]])) {
      nodes[[i]]@attrs$height <- nodeheight[[nodes[[i]]@name]]
    }
    if (!is.null(nodewidth[[nodes[[i]]@name]])) {
      nodes[[i]]@attrs$width <- nodewidth[[nodes[[i]]@name]]
    }
    if (length(grep("and", nodes[[i]]@name)) > 0) {
      if (is.null(nodelabel)) {
        nodelabel <- list()
      }
      nodelabel[[nodes[[i]]@name]] <- "AND"
      nodes[[i]]@attrs$label <- ""
      nodes[[i]]@attrs$fontcolor <- andcolor
      if (is.null(bordercol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$color <- "grey"
      } else {
        nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
      }
      if (is.null(nodeshape[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$shape <- "ellipse"
      } else {
        nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
      }
      if (is.null(nodecol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$fillcolor <- "grey"
      } else {
        nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
      }
      nodes[[i]]@attrs$width <- "0.5"
      nodes[[i]]@attrs$height <- "0.5"
        if (type == 2) {
          nodes[[i]]@attrs$fontsize <- "0"
        } else {
          nodes[[i]]@attrs$fontsize <- "0"
        }
    } else {
      nodes[[i]]@attrs$fontsize <- as.character(fontsize)
      if (is.null(nodecol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$fillcolor <- "white"
      } else {
        nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
      }
      if (is.null(nodeshape[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$shape <- "ellipse"
      } else {
        nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
      }
      if (is.null(bordercol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$color <- "black"
      } else {
        nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
      }
      if (names(nodes)[i] %in% stimuli & is.null(nodeshape[[nodes[[i]]@name]])) {
        if (type == 2) {
          nodes[[i]]@attrs$shape <- "diamond"
        } else {
          nodes[[i]]@attrs$shape <- "box"
        }
      }
      if (names(nodes)[i] %in% signals & is.null(nodecol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$fillcolor <- "lightblue"
      }
      if (names(nodes)[i] %in% inhibitors & is.null(bordercol[[nodes[[i]]@name]])) {
        nodes[[i]]@attrs$color <- "red"
      }
    }
    if (!is.null(nodestates)) {
      if (sum(names(nodestates) %in% nodes[[i]]@name) == 1) {
        if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 0) {
          nodes[[i]]@attrs$fillcolor <- "white"
          nodes[[i]]@attrs$color <- "black"
        }
        if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 1) {
          nodes[[i]]@attrs$fillcolor <- "green"
          nodes[[i]]@attrs$color <- "black"
        }
      }
    }
  }
  if (length(edges) > 0) {
    for (i in 1:length(edges)) {
      edges[[i]]@attrs$fontsize <- as.character(labelsize)
      if (length(grep("and", names(edges)[i])) > 0) {
        tmp <- unlist(strsplit(names(edges)[i], "~"))
        k <- grep("and", tmp)
        inputN <- length(grep(tmp[k], edges))
        k <- as.numeric(gsub("and", "", tmp[k]))
        inputN2 <- grep(tmp[1], unlist(strsplit(dnf[grep("\\+", dnf)[k]], "\\+")))-1
        edges[[i]]@attrs$style <- lty[grep("\\+", dnf)[k]]
        edges[[i]]@attrs$label <- labels[grep("\\+", dnf)[k]]
        if (use.freq) {
          edges[[i]]@attrs$weight <- freq[grep("\\+", dnf)[k]]
          edges[[i]]@attrs$fontcolor <- "blue"
        }
        if (!is.null(edgewidth)) {
          edges[[i]]@attrs$label <- edgewidth[grep("\\+", dnf)[k]]
        }
        if (!is.null(edgestyle)) {
          if (!is.na(edgestyle[grep("\\+", dnf)[k]])) {
            edges[[i]]@attrs$style <- edgestyle[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
          }
        }
        if (!is.null(edgelabel)) {
          if (!is.na(edgelabel[grep("\\+", dnf)[k]])) {
            edges[[i]]@attrs$label <- edgelabel[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
          }
        }
        if (length(grep("!", names(edgesneg)[i])) > 0) {
          edges[[i]]@attrs$arrowhead <- "tee"
          edges[[i]]@attrs$color <- "red"
          if (!is.null(edgecol)) {
            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
              edges[[i]]@attrs$color <- edgecol[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
            }
          }
          if (!is.null(edgehead)) {
            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
              edges[[i]]@attrs$arrowhead <- edgehead[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
            }
          }
        } else {
          edges[[i]]@attrs$arrowhead <- "open"
          edges[[i]]@attrs$color <- "green"
          if (gsub("and.*", "and", tmp[1]) %in% "and") {
            if (!is.null(edgecol)) {
              if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                edges[[i]]@attrs$color <- edgecol[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k])])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k])])) > 1)]
              }
            }
            if (!is.null(edgehead)) {
              if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                edges[[i]]@attrs$arrowhead <- edgehead[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k])])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k])])) > 1)]
              }
            }
          } else {
            if (!is.null(edgecol)) {
              if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                edges[[i]]@attrs$color <- edgecol[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
              }
            }
            if (!is.null(edgehead)) {
              if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                edges[[i]]@attrs$arrowhead <- edgehead[grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
              }
            }
          }
        }
      } else {
        tmp <- unlist(strsplit(names(edges)[i], "~"))
        if (length(grep("!", names(edgesneg)[i])) > 0) {
          edges[[i]]@attrs$style <- lty[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          edges[[i]]@attrs$label <- labels[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          if (use.freq) {
            edges[[i]]@attrs$weight <- freq[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
            edges[[i]]@attrs$fontcolor <- "blue"
          }
          if (!is.null(edgewidth)) {
            edges[[i]]@attrs$label <- edgewidth[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          }
          if (!is.null(edgestyle)) {
            if (!is.na(edgestyle[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$style <- edgestyle[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          if (!is.null(edgelabel)) {
            if (!is.na(edgelabel[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$label <- edgelabel[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          edges[[i]]@attrs$arrowhead <- "tee"
          edges[[i]]@attrs$color <- "red"
          if (!is.null(edgecol)) {
            if (!is.na(edgecol[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$color <- edgecol[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          if (!is.null(edgehead)) {
            if (!is.na(edgehead[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$arrowhead <- edgehead[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
        } else {
          edges[[i]]@attrs$style <- lty[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          edges[[i]]@attrs$label <- labels[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          if (use.freq) {
            edges[[i]]@attrs$weight <- freq[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
            edges[[i]]@attrs$fontcolor <- "blue"
          }
          if (!is.null(edgewidth)) {
            edges[[i]]@attrs$label <- edgewidth[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
          }
          if (!is.null(edgestyle)) {
            if (!is.na(edgestyle[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$style <- edgestyle[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          if (!is.null(edgelabel)) {
            if (!is.na(edgelabel[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$label <- edgelabel[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          edges[[i]]@attrs$arrowhead <- "open"
          edges[[i]]@attrs$color <- "green"
          if (!is.null(edgecol)) {
            if (!is.na(edgecol[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$color <- edgecol[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
          if (!is.null(edgehead)) {
            if (!is.na(edgehead[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
              edges[[i]]@attrs$arrowhead <- edgehead[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
            }
          }
        }
      }
    }
  }
  if (type == 1) {
  
    g2 <- agopen(name="boolean", nodes=nodes, recipEdges = "distinct", edges=edges, edgeMode="undirected", attrs=list(edge = list(), graph = list(lwd = lwd, rankdir = rankdir), node=list(lwd = lwd, fixedsize=FALSE)))

    plot(g2, "dot", lwd = lwd, ...)

  } else {

    arrowheads <- character()
    arrowcolors <- character()
    arrowlabels <- character()
    arrowlwd <- character()
    arrowlty <- character()
    arrowfontsize <- character()
    if (length(edges) > 0) {
      for (i in 1:length(edges)) {
        if (length(edges[[i]]@attrs$style) == 0) { edges[[i]]@attrs$style <- "solid" }
        arrowlty <- c(arrowlty, edges[[i]]@attrs$style)
        arrowheads <- c(arrowheads, edges[[i]]@attrs$arrowhead)
        arrowcolors <- c(arrowcolors, edges[[i]]@attrs$color)
        arrowfontsize <- c(arrowfontsize, edges[[i]]@attrs$fontsize)
        arrowlwd <- c(arrowlwd, edges[[i]]@attrs$weight)
        arrowlabels <- c(arrowlabels, edges[[i]]@attrs$label)
      }
    }
    arrowlwd <- as.numeric(arrowlwd)
    
    graph.trans <- NULL
    and.count <- 0
    for (i in dnf) {
      if (length(grep("\\+", i)) > 0) {
        and.count <- and.count + 1
        output <- unlist(strsplit(i, "="))
        input <- unlist(strsplit(output[1], "\\+"))
        output <- output[2]
        for (i in input) {
          graph.trans <- c(graph.trans, paste(gsub("!", "", i), "~", paste("and", and.count, sep = ""), sep = ""))
        }
        graph.trans <- c(graph.trans, paste(paste("and", and.count, sep = ""), "~", output, sep = ""))
      } else {
        put <- unlist(strsplit(i, "="))
        graph.trans <- c(graph.trans, paste(gsub("!", "", put[1]), "~", put[2], sep = ""))
      }
    }

    if (length(edgecol) == length(arrowcolors)) {
      edgecol <- edgecol[order(match(graph.trans, names(edges)))]
      arrowcolors <- edgecol
    }
    
    nodeshapes <- character()
    nodecolors <- character()
    nodeheight <- character()
    nodewidth <- character()
    nodecolor <- character()
    
    for (i in 1:length(nodes)) {
      nodeshapes <- c(nodeshapes, nodes[[i]]@attrs$shape)
      nodecolors <- c(nodecolors, nodes[[i]]@attrs$fillcolor)
      nodeheight <- c(nodeheight, nodes[[i]]@attrs$height)
      nodewidth <- c(nodewidth, nodes[[i]]@attrs$width)
      nodecolor <- c(nodecolor, nodes[[i]]@attrs$color)
    }

    nodeheight[which(nodeheight == "0.4")] <- "0.2"

    if (is.null(lty) & is.null(edgestyle)) {
      arrowlty <- rep("solid", length(edges))
    }
    
    if (use.freq) {
      if (is.null(lty)) {
        if (edgestyle) {
          arrowlty[which(as.numeric(arrowlabels) < 66)] <- "dashed"
          arrowlty[which(as.numeric(arrowlabels) < 33)] <- "dotted"
        }
      }
      arrowlwd <- arrowlwd - min(arrowlwd)
      arrowlwd <- as.character((arrowlwd/max(arrowlwd)+0.1)*2*edgelwd)
    } else {
      arrowlwd <- rep(edgelwd, length(edges))
    }

    if (is.null(edgewidth) & is.null(edgelabel)) {
      arrowlabels <- rep("", length(edges))
    }
    
    if (length(arrowlty) == 0) {
      arrowlty <- rep("solid", length(edges))
    }
    if (length(arrowlwd) == 0) {
      arrowlwd <- rep(lwd, length(edges))
    }
    
    names(arrowfontsize) <- names(arrowheads) <- names(arrowcolors) <- names(arrowlwd) <- names(arrowlty) <- names(arrowlabels) <- names(edges)

    names(nodecolor) <- names(nodewidth) <- names(nodeheight) <- names(nodeshapes) <- names(nodecolors) <- names(nodes)

    if (length(unique(names(edges))) < length(names(edges))) {
      for (i in names(edges)[-which(duplicated(names(edges)) == TRUE)]) {
        getpos <- grep(paste("^", i, "$", sep = ""), names(edges))
        if (length(getpos) > 1) {
          if (use.freq) {
            if (arrowheads[getpos[1]] %in% "tee") {
              arrowlabels[getpos[1]] <- paste(paste(c("-", "+"), arrowlabels[getpos], sep = ""), collapse = "\n")
            } else {
              arrowlabels[getpos[1]] <- paste(paste(c("+", "-"), arrowlabels[getpos], sep = ""), collapse = "\n")
            }
          } else {
            if (is.null(edgelabel)) {
              arrowlabels[getpos[1]] <- ""
            }
          }
          arrowheads[getpos[1]] <- "odiamond"
          if (is.null(edgecol)) {
            arrowcolors[getpos[1]] <- "black"
          } else {
            if (is.na(edgecol[getpos[1]])) {
              arrowcolors[getpos[1]] <- "black"
            }
          }
          arrowlwd[getpos[1]] <- as.character(mean(as.numeric(arrowlwd[getpos])))
        }
      }
    }
    if (length(labels) == length(arrowlabels) & is.null(edgelabel)) {
      arrowlabels[!is.na(labels)] <- labels[!is.na(labels)]
    }
    if (length(edgecol) == 1) {
      arrowcolors <- rep(edgecol, length(arrowcolors))
      names(arrowcolors) <- names(arrowlabels)
    }
    
    if (legend == 1 | legend == 3) {
      if (dolegend) {
        start <- 1
        g@nodes <- c("LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
        g@edgeL <- list()
        g@edgeData@data <- list()
      } else {
        start <- length(g@nodes) + 1
        g@nodes <- c(g@nodes, "LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
      }
      g@edgeL[["LEGEND:"]] <- list()
      g@edgeL[["STIMULUS"]] <- list()
      g@edgeL[["INHIBITOR"]] <- list()
      g@edgeL[["SIGNAL"]] <- list()
      g@edgeL[["NOTHING"]] <- list()
      g@edgeL[["active"]] <- list()
      g@edgeL[["inactive"]] <- list()
      g@edgeL[["LEGEND:"]][["edges"]] <- as.integer(start+1)
      g@edgeL[["STIMULUS"]][["edges"]] <- as.integer(start+2)
      g@edgeL[["INHIBITOR"]][["edges"]] <- as.integer(start+3)
      g@edgeL[["SIGNAL"]][["edges"]] <- as.integer(start+4)
      g@edgeL[["NOTHING"]][["edges"]] <- as.integer(start+5)
      g@edgeL[["active"]][["edges"]] <- c(as.integer(start+6), as.integer(start+4))
      g@edgeL[["inactive"]][["edges"]] <- as.integer(start+5)
      g@edgeData@data[["LEGEND:|STIMULUS"]] <- list()
      g@edgeData@data[["STIMULUS|INHIBITOR"]] <- list()
      g@edgeData@data[["INHIBITOR|SIGNAL"]] <- list()
      g@edgeData@data[["SIGNAL|NOTHING"]] <- list()
      g@edgeData@data[["NOTHING|active"]] <- list()
      g@edgeData@data[["active|inactive"]] <- list()
      g@edgeData@data[["active|NOTHING"]] <- list()
      g@edgeData@data[["inactive|active"]] <- list()
      g@edgeData@data[["LEGEND:|STIMULUS"]][["weight"]] <- 1
      g@edgeData@data[["STIMULUS|INHIBITOR"]][["weight"]] <- 1
      g@edgeData@data[["INHIBITOR|SIGNAL"]][["weight"]] <- 1
      g@edgeData@data[["SIGNAL|NOTHING"]][["weight"]] <- 1
      g@edgeData@data[["NOTHING|active"]][["weight"]] <- 1
      g@edgeData@data[["active|inactive"]][["weight"]] <- 1
      g@edgeData@data[["active|NOTHING"]][["weight"]] <- 1
      g@edgeData@data[["inactive|active"]][["weight"]] <- 1
      arrowheads <- c(arrowheads, "LEGEND:~STIMULUS" = "none", "STIMULUS~INHIBITOR" = "open", "INHIBITOR~SIGNAL" = "tee", "SIGNAL~NOTHING" = "odiamond", "NOTHING~active" = "open", "active~inactive" = "tee", "active~NOTHING" = "tee", "inactive~active" = "open")
      arrowcolors <- c(arrowcolors, "LEGEND:~STIMULUS" = "transparent", "STIMULUS~INHIBITOR" = "green", "INHIBITOR~SIGNAL" = "red", "SIGNAL~NOTHING" = "black", "NOTHING~active" = "green", "active~inactive" = "red", "active~NOTHING" = "red", "inactive~active" = "green")
      arrowlabels <- c(arrowlabels, "LEGEND:~STIMULUS" = "", "STIMULUS~INHIBITOR" = "positive", "INHIBITOR~SIGNAL" = "negative", "SIGNAL~NOTHING" = "ambiguous\npositive\nnegative", "NOTHING~active" = "bidirectional\ndifferent", "active~inactive" = "bidirectional\ndifferent", "active~NOTHING" = "", "inactive~active" = "")
      nodecolors <- c(nodecolors, "LEGEND:" = "white", "STIMULUS" = "white", "INHIBITOR" = "white", "SIGNAL" = "lightblue", "NOTHING" = "white", "active" = "green", "inactive" = "white")
      nodeheight <- c(nodeheight, "LEGEND:" = 0, "STIMULUS" = as.character(max(nodeheight)), "INHIBITOR" = as.character(max(nodeheight)), "SIGNAL" = as.character(max(nodeheight)), "NOTHING" = as.character(max(nodeheight)), "active" = as.character(max(nodeheight)), "inactive" = as.character(max(nodeheight)))
      nodewidth <- c(nodewidth, "LEGEND:" = as.character(max(nodewidth)), "STIMULUS" = as.character(max(nodewidth)), "INHIBITOR" = as.character(max(nodewidth)), "SIGNAL" = as.character(max(nodewidth)), "NOTHING" = as.character(max(nodewidth)), "active" = as.character(max(nodewidth)), "inactive" = as.character(max(nodewidth)))
      if (type == 2) {
        nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "diamond", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
      } else {
        nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "box", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
      }
      nodecolor <- c(nodecolor, "LEGEND:" = "white", "STIMULUS" = "black", "INHIBITOR" = "red", "SIGNAL" = "white", "NOTHING" = "black", "active" = "black", "inactive" = "black")
      dnf <- c(dnf, "NOTHING=active", "!active=NOTHING", "!active=inactive", "inactive=active")
    }
    nodelabels <- names(nodecolor)
    names(nodelabels) <- nodelabels
    for (i in 1:length(nodelabel)) {
      nodelabels[which(names(nodelabels) %in% names(nodelabel)[i])] <- nodelabel[i]
    }
    nodefontsizes <- NULL
    if (!is.null(nodefontsize)) {
      nodefontsizes <- rep(14, length(nodelabels))
      names(nodefontsizes) <- names(nodelabels)
      for (i in 1:length(nodefontsize)) {
        nodefontsizes[which(names(nodefontsizes) %in% names(nodefontsize)[i])] <- nodefontsize[[i]]
      }
    }
    g <- layoutGraph(g, edgeAttrs = list(arrowhead = arrowheads, color = arrowcolors, label = arrowlabels), nodeAttrs = list(labels = nodelabels, color = nodecolor, height = nodeheight, width = nodewidth, shape = nodeshapes, fillcolor = nodecolors), layoutType=layout)
    graph.par(list(graph=list(main = main, sub = sub, cex.main = cex.main, cex.sub = cex.sub, col.sub = col.sub), edges=list(textCol = labelcol, lwd = edgelwd, fontsize = labelsize), nodes=list(lwd = lwd, fontsize = fontsize, cex = cex)))
    edgeRenderInfo(g) <- list(lty = arrowlty, lwd = arrowlwd, label = arrowlabels)
    if (length(edges) > 0) {
      for (i in names(g@renderInfo@edges$direction)) {
        if (length(grep("and", i)) == 0 & g@renderInfo@edges$direction[[i]] == "both") {
          pos <- which(names(g@renderInfo@edges$arrowhead) %in% i)
          input <- unlist(strsplit(i, "~"))
          output <- input[2]
          input <- input[1]
          if (paste("!", input, "=", output, sep = "") %in% dnf) {
            g@renderInfo@edges$arrowhead[pos] <- "tee"
          } else {
            g@renderInfo@edges$arrowhead[pos] <- "open"
          }
          if (paste("!", output, "=", input, sep = "") %in% dnf) {
            g@renderInfo@edges$arrowtail[pos] <- "tee"
          } else {
            g@renderInfo@edges$arrowtail[pos] <- "open"
          }
          if (paste("!", output, "=", input, sep = "") %in% dnf & paste("", output, "=", input, sep = "") %in% dnf) {
            g@renderInfo@edges$arrowtail[pos] <- "odiamond"
          }
          if (paste("!", input, "=", output, sep = "") %in% dnf & paste("", input, "=", output, sep = "") %in% dnf) {
            g@renderInfo@edges$arrowhead[pos] <- "odiamond"
          }
          if (is.null(edgecol)) {
            if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
              g@renderInfo@edges$col[pos] <- "green"
            }
            if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
              g@renderInfo@edges$col[pos] <- "red"
            }
            if (g@renderInfo@edges$arrowtail[pos] == "odiamond" & g@renderInfo@edges$arrowhead[pos] == "odiamond") {
              g@renderInfo@edges$col[pos] <- "black"
            }
            if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
              g@renderInfo@edges$col[pos] <- "brown"
            }
          } else {
            if (is.null(edgecol)) { # is.na(edgecol[pos])
              if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
                g@renderInfo@edges$col[pos] <- "green"
              }
              if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
                g@renderInfo@edges$col[pos] <- "red"
              }
              if (g@renderInfo@edges$arrowtail[pos] == "odiamond" & g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                g@renderInfo@edges$col[pos] <- "black"
              }
              if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
                g@renderInfo@edges$col[pos] <- "brown"
              }
            }
          }
        }
      }
    }
    #g@renderInfo@nodes$labelX[grep("and", names(g@renderInfo@nodes$labelX))] <- -1000
    #g@renderInfo@nodes$labelY[grep("and", names(g@renderInfo@nodes$labelY))] <- -1000
    renderGraph(g, lwd = lwd, ...)
  }

  if (legend == 2 | legend == 3) {
    legend(x = x, y = y, legend = c("signals are blue", "stimuli are diamonds/boxes", "inhibitors have a red border", "positive regulation is green ->", "negative regulation is red -|", "ambiguous regulation is black -o"), fill = c("lightblue", "white", "red", "green", "red", "black"), col = c("lightblue", "white", "red", "green", "red", "black"), yjust = yjust, xjust = xjust)
  }
  return(g)
  
}


plotMovie <- function(graph, stimuli, save = NULL, sleep = 1, ...) {
  count <- 1
  hierarchy <- getHierarchy(graph)
  nodes <- unlist(hierarchy)
  if (!is.null(save)) {
    pdf(paste("temp", count, ".pdf", sep = ""), ...)
  }
  plotDnf(graph, simulate = list(stimuli = stimuli, inhibitors = nodes[-which(nodes %in% stimuli)]), ...)
  count <- count + 1
  if (!is.null(save)) {
    dev.off()
  }
  Sys.sleep(sleep)
  nodes.done <- NULL
  for (i in 2:length(hierarchy)) {
    nodes.done <- c(nodes.done, hierarchy[[i]])
    if (!is.null(save)) {
      pdf(paste("temp", count, ".pdf", sep = ""), ...)
    }
    plotDnf(graph, simulate = list(stimuli = stimuli, inhibitors = hierarchy[[i]]), ...)
  count <- count + 1
    if (!is.null(save)) {
      dev.off()
    }
    Sys.sleep(sleep)
  }
  if (!is.null(save)) {
    pdf(paste("temp", count, ".pdf", sep = ""), ...)
  }
  plotDnf(graph, simulate = list(stimuli = stimuli), ...)
  count <- count + 1
  if (!is.null(save)) {
    dev.off()
  }
}

plotMovieAwsome <- function(graph, stimuli, file = "temp.gif", ...) {
  hierarchy <- getHierarchy(graph)
  nodes <- unlist(hierarchy)
  require(animation)
  saveGIF({
    plotDnf(graph, simulate = list(stimuli = stimuli, inhibitors = nodes[-which(nodes %in% stimuli)]), ...)
    nodes.done <- NULL
    for (i in 2:length(hierarchy)) {
      nodes.done <- c(nodes.done, hierarchy[[i]])
      plotDnf(graph, simulate = list(stimuli = stimuli, inhibitors = hierarchy[[i]]), ...)
    }
    plotDnf(graph, simulate = list(stimuli = stimuli), ...)
  }, movie.name = file)
}

plotDnfOld <- function(dnf, bString = NULL, CNOlist = NULL, showall = FALSE, ...) {
  dnf2 <- dnf[grep("=", dnf)]
  model <- list()
  vertices <- unlist(strsplit(dnf, "="))
  vertices <- unique(gsub("!", "", unlist(strsplit(vertices, "\\+"))))
  dnf <- dnf2
  model$reacID <- dnf
  model$interMat <- matrix(0, length(vertices), length(dnf))
  model$notMat <- matrix(0, length(vertices), length(dnf))
  colnames(model$interMat) <- dnf
  colnames(model$notMat) <- dnf
  rownames(model$interMat) <- vertices
  rownames(model$notMat) <- vertices
  model$newANDs <- as.list(dnf[grep("\\+", dnf)])
  model$speciesCompressed <- NULL
  model$SplitANDs$initialReac <- c("split1", "split2")
  for (i in dnf) {
    tmp <- unlist(strsplit(i, "="))
    tmp2 <- unlist(strsplit(tmp[1], "\\+"))
    model$interMat[which(rownames(model$interMat) %in% tmp[2]), which(colnames(model$interMat) %in% i)] <- 1
    for (j in tmp2) {
      model$interMat[which(rownames(model$interMat) %in% gsub("!", "", j)), which(colnames(model$interMat) %in% i)] <- -1
      if (!(gsub("!", "", j) %in% j)) {
        model$notMat[which(rownames(model$notMat) %in% gsub("!", "", j)), which(colnames(model$notMat) %in% i)] <- 1
      }
    }
  }
  if (is.null(CNOlist) | !(showall)) {
    model$namesSpecies <- vertices
  } else {
    model$namesSpecies <- unique(c(vertices, colnames(CNOlist@signals[[1]])))
  }
  plotModel(rep(1, length(model$reacID)), model = model, CNOlist = CNOlist, bString = bString, indexIntegr = NULL, ...)
}
