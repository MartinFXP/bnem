residualErrorNem <- function(CNOlist, NEMlist, model, approach = "fc", sizeFac = 1e-4, NAFac = 1, parameters = c(1, 0.2, 1, 2)) {
  
  NEMlist <- computeFc(CNOlist, NEMlist$exprs, type = "data", parameters = parameters)
  
  dataCont <- NEMlist$fc

  dataCont[which(abs(dataCont) <= parameters[2])] <- 0

  dataCont[which(dataCont >= parameters[2])] <- 1

  dataCont[which(dataCont <= -parameters[2])] <- -1
  
  dataDisc <- round(dataCont)

  characters <- numeric()

  for (i in 1:nrow(dataNorm)) {

    characters <- c(characters, paste(as.character(dataNorm[i, ]), collapse=""))

  }

  uniques <- unique(characters)

  clusters <- numeric()

  for (i in 1:length(uniques)) {
    clusters[i] <- sum(characters %in% uniques[i])
  }

  bigClusts <- uniques[order(clusters, decreasing = TRUE)[1:ncol(CNOlist@signals[[1]])]]

  Sgenes <- numeric()
  
  ## for (i in bigClusts) {
  ##   Sgenes <- rbind(Sgenes, 
  
}

