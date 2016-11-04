kmeansNorm <-
function(x, k = 2) {
  if (!is.matrix(x)) {
    x <- t(as.matrix(x))
    if (nrow(x) != 1) {
      x <- t(x)
    }
  }
  require(cluster)
  y <- x
  for (i in 1:nrow(x)) {

    if (sd(x[i, ]) == 0) { next() }
    
    cat('\r', paste(round(i/nrow(x)*100), "%", sep = ""))
    flush.console()
        
    x.clust <- kmeans(x[i, ], k)
    x.dist <- dist(x[i, ])
    x.sil <- silhouette(x.clust$cluster, x.dist)
    x.norm <- x.sil[, 3]
    if (sum(x[i, which(x.clust$cluster == 1)]) < sum(x[i, which(x.clust$cluster == 2)])) {
      x.norm[which(x.clust$cluster == 1)] <- x.norm[which(x.clust$cluster == 1)]*(-1)
    } else {
      x.norm[which(x.clust$cluster == 2)] <- x.norm[which(x.clust$cluster == 2)]*(-1)
    }
    x.norm <- x.norm - min(x.norm)
    x.norm <- x.norm/max(x.norm)
    y[i, ] <- x.norm
  }
  return(y)
}
