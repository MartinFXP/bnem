simpleNorm <- function(x) {
  if (is.matrix(x) == FALSE) {
    genes.min <- min(x)
    x <- (x - genes.min)
    genes.max <- max(x)
    x <- x/genes.max
  } else {
    genes.min <- apply(x, 1, min)
    x <- (x - genes.min)
    genes.max <- apply(x, 1, max)
    x <- x/genes.max
  }
  x[is.na(x)] <- 0
  return(x)
}

kmeansNorm <- function(x) {
  require(cluster)
  y <- x
  for (i in 1:nrow(x)) {
    x.clust <- kmeans(x[i, ], 2)
    x.dist <- dist(x[i, ])
    x.sil <- silhouette(x.clust$cluster, x.dist)
    x.norm <- x.sil
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

disc <- function(x, beta = 0.5) { # simple discretization
  x.disc <- x
  x.disc[which(abs(x.disc) <= beta)] <- 0
  x.disc[which(x.disc > beta)] <- 1
  x.disc[which(x.disc < -beta)] <- -1
  return(x.disc)
}

pamNorm <- function(x, method = "euclidean") { # 2pam clustering and silhoutte normalization
  require(cluster)
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

expNorm <- function(x,  stimuli = NULL, inhibitors = NULL, batches, runs, cutoff) {
  design <- makeDesign(x, stimuli, inhibitors, c(batches, runs))
  # for every gene find all the batches in all the runs that have a significant change and normalize (pam or simple)
  # normalize batches that do not have a significant change according to a level cutoff or to control level
  normedX <- x*0
  for (run in runs) {
    for (batch in batches) {
      targetRows <- intersect(which(design[, run] == 1), which(design[, batch] == 1))
      if (length(targetRows) == 0) { next() }
      for (i in 1:nrow(normedX)) {
        if (max(x[i, targetRows]) - min(x[i, targetRows]) >= cutoff) {
          normedX[i, targetRows] <- simpleNorm(x[i, targetRows])
        }
      }
    }
  }
  cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors), collapse = "|"), colnames(design))], 1, sum)
  grepCtrl <- which(cuesSum == 0)
  for (run in runs) {
    for (batch in batches) {
      targetRows <- intersect(which(design[, run] == 1), which(design[, batch] == 1))
      if (length(targetRows) == 0) { next() }
      for (i in 1:nrow(normedX)) {
        if (max(x[i, targetRows]) - min(x[i, targetRows]) < cutoff) {
          normedX[i, targetRows] <- median(normedX[i, grepCtrl])
        }
      }
    }
  }
  return(normedX)
}
