heatmapOP <-
function(x, col = "RdYlGn", coln = 11, bordercol = "grey", borderwidth = 0.1, breaks = "sym", main = "heatmap by Oscar Perpiñán", sub = "http://oscarperpinan.github.io/rastervis/; http://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps", dendrogram = "none", colorkey = list(space = "right"), Colv = TRUE, Rowv = TRUE, xrot = 90, yrot = 0, shrink = c(1,1), cexCol = 1, cexRow = 1, cexMain = 1, cexSub = 1, colSideColors = NULL, aspect = "fill", contour = FALSE, useRaster = FALSE, xlab = NULL, ylab = NULL, colSideColorsPos = "top", clust = NULL, clusterx = NULL, ...) { ## if clust == "p","s","k" uses correlation to cluster with clust as method
  ## for info on the commands checkout ?levelplot and ?xyplot in library(lattice)
  ## heatmap © Oscar Perpiñán @ http://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps; http://oscarperpinan.github.io/rastervis/
  ## x - matrix
  ## col - color from brewer.pal.info
  ## bordercol - normal colorname or number
  ## borderwidth
  ## breaks - either number or numeric vector for the breaks
  ## main
  ## sub
  ## dendrogram - "both", "col" or "row"
  ## colorkey - list(space = "right") or "top", "bottom", "left"
  ## Colv - true or false for clustering of columns
  ## Rowv - true or false for clustering of rows
  ## xrot - rotation of the column names
  ## yrot - rotation of the rownames
  ## shrink - numeric vector for size of the fields (in percent) of the small to the large numbers
  ## cexCol,cexRow,cexMain,cexSub - fontsize
  ## aspect = "fill" for streched squares, or "iso" for quadratic ones
  if (sum(is.na(x)) > 0) {
    print("NAs detected; set to 0")
    x[is.na(x)] <- 0
  }

  if (max(x) == min(x)) {
    x <- matrix(rnorm(100), 10, 10)
    main <- "max value equals min value"
    sub <- "random matrix plotted"
  }

  dd.col <- NULL

  if (dendrogram == "row" | dendrogram == "both" & !is.null(colorkey)) {
    colorkey = list(space = "left")
  }

  if (is.null(breaks)) {
    breaks <- seq(min(x),max(x),(max(x) - min(x))/45)
  }
  if ("sym" %in% breaks) {
    breaks <- max(abs(x))
    breaks2 <- breaks/45
    breaks <- seq(-breaks,breaks,breaks2)
  }
  if (length(breaks) == 1) {
    at <- seq(min(x), max(x), (max(x)-min(x))/breaks)
  } else {
    at <- breaks
    x[x < breaks[1]] <- breaks[1]
    x[x > breaks[length(breaks)]] <- breaks[length(breaks)]
  }
  
  require(lattice)
  require(latticeExtra)
  
  if (Colv) {
    if (is.null(clust)) {
      if (is.null(clusterx)) {
        dd.col <- as.dendrogram(hclust(dist(t(x))))
      } else {
        dd.col <- as.dendrogram(hclust(dist(t(clusterx))))
      }
    } else {
      if (is.null(clusterx)) {
        cor.tmp <- cor(x, method = clust)
      } else {
        cor.tmp <- cor(clusterx, method = clust)
      }
      cor.tmp[which(is.na(cor.tmp) == TRUE)] <- 0
      dd.col <- as.dendrogram(hclust(as.dist(1-abs(cor.tmp))))
    }
    col.ord <- order.dendrogram(dd.col)
  } else {
    if (dendrogram %in% "both") {
      dendrogram <- "row"
    }
    if (dendrogram %in% "col") {
      dendrogram <- "none"
    }
    col.ord <- 1:ncol(x)
    legend = NULL
  }

  if (Rowv) {
    if (is.null(clust)) {
      if (is.null(clusterx)) {
        dd.row <- as.dendrogram(hclust(dist(x)))
      } else {
        dd.row <- as.dendrogram(hclust(dist(clusterx)))
      }
    } else {
      if (is.null(clusterx)) {
        cor.tmp <- cor(t(x), method = clust)
      } else {
        cor.tmp <- cor(t(clusterx), method = clust)
      }
      cor.tmp[which(is.na(cor.tmp) == TRUE)] <- 0
      dd.row <- as.dendrogram(hclust(as.dist(1-abs(cor.tmp))))
    }
    row.ord <- rev(order.dendrogram(dd.row))
  } else {
    if (dendrogram %in% "both") {
      dendrogram <- "col"
    }
    if (dendrogram %in% "row") {
      dendrogram <- "none"
    }
    row.ord <- 1:nrow(x)
    legend = NULL
  }

  ## is the following correct? and do I have to delete this before clsutering? without this, the clustering uses, what the eye sees...
  ## if (length(breaks) > 1) {
  ##   x[x < breaks[1]] <- breaks[1]
  ##   x[x > breaks[length(breaks)]] <- breaks[length(breaks)]
  ## }
  
  ##Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds
  ##YlGn YlGnBu YlOrBr YlOrRd
  ##All the sequential palettes are available in variations from 3 different values up to 9 different values.
  ##The diverging palettes are
  ##BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
  ## brewer.pal.info for all of them

  add <- list(rect = list(col = "transparent", fill = colSideColors[sort(col.ord)]))

  myTheme <- custom.theme(region=brewer.pal(n=coln, col))
  if (dendrogram != "none") {
    if (dendrogram == "both") {
      if (colSideColorsPos %in% "bottom") {
        if (!is.null(colSideColors)) {
          size <- 2
        } else {
          size <- 10
        }
        legend <- list(bottom =
                       list(fun = dendrogramGrob,
                            args =
                            list(x = dd.col, ord = col.ord,
                                 side = "top",
                                 size = size,
                                 add = add,
                                 type = "rectangle")),
                       right =
                       list(fun = dendrogramGrob,
                            args =
                            list(x = dd.row,
                                 side = "right", size = 10)))
      }
      if (colSideColorsPos %in% "top") {
        legend <- list(top =
                       list(fun = dendrogramGrob,
                            args =
                            list(x = dd.col, ord = col.ord,
                                 side = "top",
                                 size = 10,
                                 add = add,
                                 type = "rectangle")),
                       right =
                       list(fun = dendrogramGrob,
                            args =
                            list(x = dd.row,
                                 side = "right", size = 10)))
      }
    }
    if (dendrogram == "row") {
      if (!is.null(colSideColors)) {
        if (is.null(dd.col)) {
          col.ord <- 1:length(colSideColors)
        }
        if (is.null(clusterx)) {
          dd.col <- as.dendrogram(hclust(dist(t(x))*0))
        } else {
          dd.col <- as.dendrogram(hclust(dist(t(clusterx))*0))
        }
        if (colSideColorsPos %in% "bottom") {
          legend <- list(right =
                         list(fun = dendrogramGrob,
                              args =
                              list(x = dd.row,
                                   side = "right",
                                   size = 10, size.add= 0.5)),
                         bottom =
                         list(fun = dendrogramGrob,
                              args = list(x = dd.col, ord = col.ord, 
                                side = "top", size = 1, size.add= 1,
                                add = add,
                                type = "rectangle")))
        }
        if (colSideColorsPos %in% "top") {
          legend <- list(right =
                         list(fun = dendrogramGrob,
                              args =
                              list(x = dd.row,
                                   side = "right",
                                   size = 10, size.add= 0.5)),
                         top =
                         list(fun = dendrogramGrob,
                              args = list(x = dd.col, ord = col.ord, 
                                side = "top", size = 1, size.add= 1,
                                add = add,
                                type = "rectangle")))
        }
      } else {
        legend <- list(right =
                       list(fun = dendrogramGrob,
                            args =
                            list(x = dd.row,
                                 side = "right",
                                 size = 10, size.add= 0.5)))
      }
    }
    if (dendrogram == "col") {
      if (colSideColorsPos %in% "bottom") {
        if (!is.null(colSideColors)) {
          size <- 2
        } else {
          size <- 10
        }
        legend <- list(bottom =
                       list(fun = dendrogramGrob,
                            args = list(x = dd.col, ord = col.ord, 
                              side = "top", size = size, size.add= 0.5,
                              add = add,
                              type = "rectangle")))
      }
      if (colSideColorsPos %in% "top") {
        legend <- list(top =
                       list(fun = dendrogramGrob,
                            args = list(x = dd.col, ord = col.ord, 
                              side = "top", size = 10, size.add= 0.5,
                              add = add,
                              type = "rectangle")))
      }
    }
  } else {
    if (!is.null(colSideColors)) {
      if (is.null(dd.col)) {
        col.ord <- 1:length(colSideColors)
      }
      if (is.null(clusterx)) {
        dd.col <- as.dendrogram(hclust(dist(t(x))*0))
      } else {
        dd.col <- as.dendrogram(hclust(dist(t(clusterx))*0))
      }
      if (colSideColorsPos %in% "bottom") {
        legend <- list(bottom =
                       list(fun = dendrogramGrob,
                            args = list(x = dd.col, ord = col.ord, 
                              side = "top", size = 1, size.add= 1,
                              add = add,
                              type = "rectangle"))
                       )
      }
      if (colSideColorsPos %in% "top") {
        legend <- list(top =
                       list(fun = dendrogramGrob,
                            args = list(x = dd.col, ord = col.ord, 
                              side = "top", size = 1, size.add= 1,
                              add = add,
                              type = "rectangle"))
                       )
      }
    } else {
      legend <- NULL
    }
  }
  
  d <- t(x[row.ord, col.ord])

  d <- d[, ncol(d):1]

  if (contour) {
    region <- TRUE
    col.regions <- terrain.colors(100)
  }

  ##  print(p, position=c(0,ypct-0.05,1,1), more=TRUE)
    ##  print(p2, position=c(0,0,1,ypct+0.05))

    if (is.null(rownames(x))) {
        ytck <- list(cex = 0, rot = 0, at = NULL)
    } else {
        ytck <- list(cex = cexRow, rot = yrot)
    }

    if (is.null(colnames(x))) {
        xtck <- list(cex = 0, rot = 0, at = NULL)
    } else {
        xtck <- list(cex = cexCol, rot = xrot)
    }

  levelplot(d,
            main = list(label = main, cex = cexMain), 
            sub = list(label = sub, cex = cexSub),
            aspect = aspect, xlab=xlab, ylab=ylab,
            scales = list(x = xtck, y = ytck, tck = c(1,0)),
            par.settings=myTheme,
            border=bordercol, border.lwd=borderwidth,
            shrink=shrink,
            legend = legend,
            at = at,
            colorkey = colorkey,
            contour = contour)
}
