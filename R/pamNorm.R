#' @noRd
#' @import cluster
pamNorm <-
function(x, method = "euclidean") { # 2pam clustering and silhoutte normalization
  if (is.matrix(x)) {
    for (i in 1:nrow(x)) {
      clust <- pam(x[i, ], 2)
      clust.nums <- clust$clustering
      dist <- dist(x[i, ], method = method)
      sil <- silhouette(clust.nums, dist=dist)
      max <- which(x[i, ] == max(x[i, ]))
      cluster.one <- sil[max, 1]
      cluster.zero <- c(1,2)[-cluster.one]
      min <- which(x[i, ] == min(x[i, ]))
      sil[which(sil[, 1] == cluster.zero), 3] <- sil[which(sil[, 1] == cluster.zero), 3]*(-1)
      sil[, 3] <- (sil[, 3] + 1)/2
      sil[which(x[i, ] >= x[i, which(sil[, 3] == max(sil[, 3]))[1]]), 3] <- sil[which(x[i, ] == x[i, which(sil[, 3] == max(sil[, 3]))[1]])[1], 3]
      sil[which(x[i, ] <= x[i, which(sil[, 3] == min(sil[, 3]))[1]]), 3] <- sil[which(x[i, ] == x[i, which(sil[, 3] == min(sil[, 3]))[1]])[1], 3]
      sil[, 3] <- sil[, 3] - min(sil[, 3])
      sil[, 3] <- sil[, 3]/max(sil[, 3])
      x[i, ] <- sil[, 3]
    }
  } else {
    clust <- pam(x, 2)
    clust.nums <- clust$clustering
    dist <- dist(x, method = method)
    sil <- silhouette(clust.nums, dist=dist)
    max <- which(x == max(x))
    cluster.one <- sil[max, 1]
    cluster.zero <- c(1,2)[-cluster.one]
    min <- which(x == min(x))
    sil[which(sil[, 1] == cluster.zero), 3] <- sil[which(sil[, 1] == cluster.zero), 3]*(-1)
    sil[, 3] <- (sil[, 3] + 1)/2
    sil[which(x >= x[which(sil[, 3] == max(sil[, 3]))[1]]), 3] <- sil[which(x == x[which(sil[, 3] == max(sil[, 3]))[1]])[1], 3]
    sil[which(x <= x[which(sil[, 3] == min(sil[, 3]))[1]]), 3] <- sil[which(x == x[which(sil[, 3] == min(sil[, 3]))[1]])[1], 3]
    sil[, 3] <- sil[, 3] - min(sil[, 3])
    sil[, 3] <- sil[, 3]/max(sil[, 3])
    x <- sil[, 3]
  }
  return(x)
}
