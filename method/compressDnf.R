compressDnf <- function(dnf, sgenes, iter.max = NULL) {
  nodes <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))
  unseen <- nodes[-which(nodes %in% sgenes)]
  if (is.null(iter.max)) {
    iter.max <- length(unseen)
  }
  if (length(unseen) > 0) {
    count <- 0
    while(count < iter.max) {
      for (i in unseen) {
        tmp <- NULL
        edges <- dnf[grep(paste("=", i, "$", sep = ""), dnf)]
        if (sum(dnf %in% edges) > 0) {
          dnf <- dnf[-which(dnf %in% edges)]
        }
        parents <- unique(unlist(strsplit(gsub("=.*", "", edges), "\\+")))
        edges <- dnf[grep(paste(i, "=.*", sep = ""), dnf)]
        if (sum(dnf %in% edges) > 0) {
          dnf <- dnf[-which(dnf %in% edges)]
        }
        if (length(grep("^!", edges)) > 0) {
          childrenPos <- gsub(".*=", "", edges[-grep("^!", edges)])
          childrenNeg <- gsub(".*=", "", edges[grep("^!", edges)])
        } else {
          childrenPos <- gsub(".*=", "", edges)
          childrenNeg <- NULL
        }
        for (j in parents) {
          for (k in childrenPos) {
            tmp <- c(tmp, paste(j, k, sep = "="))
          }
          for (k in childrenNeg) {
            tmp <- c(tmp, paste("!", j, k, sep = "="))
          }
        }
        tmp <- gsub("!=", "!", tmp)
        dnf <- unique(c(dnf, tmp))
      }
      count <- count + 1
    }
    dnf <- gsub("!!", "", dnf)
    if (length(grep(paste(paste(unseen, "$", sep = ""), collapse = "|"), dnf)) > 0) {
      dnf <- dnf[-grep(paste(paste(unseen, "$", sep = ""), collapse = "|"), dnf)]
    }
    if (length(grep("!=", dnf)) > 0) {
      dnf <- dnf[-grep("!=", dnf)]
    }
    dnf <- unique(dnf)
  }
  return(dnf)
}
