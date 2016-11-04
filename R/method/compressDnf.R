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
        if (length(grep(paste(i, i, sep = "="), dnf)) > 0) {
          dnf <- dnf[-grep(paste(i, i, sep = "="), dnf)]
        }
        parents <- unlist(strsplit(gsub("!|=.*", "", dnf[grep(paste("=", i, sep = ""), dnf)]), "\\+"))
        if (length(parents) == 0 & length(grep(paste("\\+", i, "=|\\+", i, "\\+|^", i, "\\+|^", i, "=", sep = ""), dnf)) > 0) {
          dnf <- dnf[-grep(paste("\\+", i, "=|\\+", i, "\\+|^", i, "\\+|^", i, "=", sep = ""), dnf)]
          next()
        }
        tmp <- NULL
        edges <- dnf[grep(paste("=", i, "$", sep = ""), dnf)]
        for (j in edges) {
          for (k in dnf[grep(paste("+", i, "+.*=|+", i, "=|^", i, "=|^", i, "+", sep = ""), dnf)]) {
            dnf <- unique(c(dnf, gsub(i, gsub("=.*", "", j), dnf[grep(paste("+", i, "+.*=|+", i, "=|^", i, "=|^", i, "+", sep = ""), dnf)])))
          }
          for (k in dnf[grep(paste("!", i, "=|!", i, "\\+", sep = ""), dnf)]) {
            tmp <- unlist(strsplit(gsub("=.*", "", k), "\\+"))
            child <- gsub(".*=", "", k)
            for (l in tmp) {
              dnf <- unique(c(dnf, paste("!", l, "=", child, sep = "")))
            }
          }
        }
        if (length(grep(paste("=", i, "$", sep = ""), dnf)) > 0) {
          dnf <- dnf[-grep(paste("=", i, "$", sep = ""), dnf)]
        }
        if (length(grep(paste("^", i, "=|^!", i, "=", sep = ""), dnf)) > 0) {
          dnf <- dnf[-grep(paste("^", i, "=|^!", i, "=", sep = ""), dnf)]
        }
        if (length(grep(paste(i, i, sep = "="), dnf)) > 0) {
          dnf <- dnf[-grep(paste(i, i, sep = "="), dnf)]
        }
        parents <- unlist(strsplit(gsub("!|=.*", "", dnf[grep(paste("=", i, sep = ""), dnf)]), "\\+"))
        if (length(parents) == 0 & length(grep(paste("\\+", i, "=|\\+", i, "\\+|^", i, "\\+|^", i, "=|\\+!", i, "=|\\+!", i, "\\+|^!", i, "\\+|^!", i, "=", sep = ""), dnf)) > 0) {
          dnf <- dnf[-grep(paste("\\+", i, "=|\\+", i, "\\+|^", i, "\\+|^", i, "=|\\+!", i, "=|\\+!", i, "\\+|^!", i, "\\+|^!", i, "=", sep = ""), dnf)]
          next()
        }
      }
      if (length(grep(paste(unseen, collapse = "|"), dnf)) == 0) {
        count <- iter.max
      } else {
        count <- count + 1
      }
    }
    dnf <- gsub("!!", "", dnf)
    if (length(grep(paste(paste("=", unseen, "$", sep = ""), collapse = "|"), dnf)) > 0) {
      dnf <- dnf[-grep(paste(paste("=", unseen, "$", sep = ""), collapse = "|"), dnf)]
    }
    if (length(grep("!=", dnf)) > 0) {
      dnf <- dnf[-grep("!=", dnf)]
    }
    if (length(grep(paste(unseen, collapse = "|"), dnf)) > 0) {
      dnf <- dnf[-grep(paste(unseen, collapse = "|"), dnf)]
    }
    dnf <- unique(dnf)
  }
  return(dnf)
}

## example:

## par(mfrow=c(1,2))

## graph <- c("TNFR1=ROS", "TNFR1=NF-kB", "TNFR1=Caspases", "TNFR1=ASK1,MEKK1", "ROS=ASK1,MEKK1", "ASK1,MEKK1=MKK7", "Caspases=Apoptosis", "MKK7=JNK", "JNK=Apoptosis")

## plotDnf(graph[-c(3,7)], bordercol = "transparent", hierarchy = hierarchy, labelcol = "red", showall = F)

## plotDnf(graph[-c(3,7)], bordercol = "transparent", hierarchy = hierarchy, labelcol = "red", showall = T)

compressDnfOld <- function(dnf, sgenes, iter.max = NULL) {
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
