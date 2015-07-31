optScore <- function(NEMlist, parameters = NULL, method = "s") {

  method <- checkMethod(method)
  
  Egene.dist <- matrix(0, nrow(NEMlist$fc), (bin2dec(rep(1, ncol(NEMlist$fc)))-1))
  Egene.distI <- matrix(0, nrow(NEMlist$fc), (bin2dec(rep(1, ncol(NEMlist$fc)))-1))
  if (any(c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "spearman", "pearson", "kendall") %in% method)) {
    if ("spearman" %in% method) {
      ## NEMlist$fc <- t(apply(NEMlist$fc, 1, rank))
      ## for (i in 1:(bin2dec(rep(1, ncol(NEMlist$fc)))-1)) {
      ##   bin <- dec2bin(i)
      ##   for (j 1:(bin2dec(rep(1, sum(bin == 1)))-1)) {
      ##     bin[which(dec2bin(j) == 1)] <- -1
      ##     Egene.dist[, i] <- cor(c(rep(0, ncol(NEMlist$fc) - length(dec2bin(i))), dec2bin(i)), t(NEMlist$fc))
      ##   }
      ## }
    }
    if ("pearson" %in% method | "kendall" %in% method) {
      for (i in 1:(bin2dec(rep(1, ncol(NEMlist$fc)))-1)) {
        Egene.dist[, i] <- cor(c(rep(0, ncol(NEMlist$fc) - length(dec2bin(i))), dec2bin(i)), t(NEMlist$fc), method = method)
      }
    }
    if (length(grep("euclidean|maximum|manhattan|canberra|binary|minkowski", method)) > 0) {
      max.dist <- dist(rbind(rep(max(abs(NEMlist$fc)), ncol(NEMlist$fc)), rep(-1, ncol(NEMlist$fc))), method = method)
      require(flexclust)
      for (i in 1:(bin2dec(rep(1, ncol(NEMlist$fc)))-1)) {
        Sgene <- c(rep(0, ncol(NEMlist$fc) - length(dec2bin(i))), dec2bin(i))
        Egene.dist[, i] <- -(max.dist - dist2(NEMlist$fc, Sgene, method = method))/max.dist
        Egene.distI[, i] <- -(max.dist - dist2(-NEMlist$fc, Sgene, method = method))/max.dist
      }
      Egene.dist <- cbind(Egene.dist, Egene.distI)
    }
  } else {
    
  }
  
  tmp <- abs(Egene.dist)
  
  tmp <- sum(apply(tmp, 1, max))/nrow(NEMlist$fc)
  
  return(1 - tmp)
}
