#' plotting the observed response scheme of an effect reporter and the expected response scheme of the regulating signalling gene
 #' @param CNOlist CNOlist object.
#' @param fc matrix of foldchanges (observed response scheme, ORS).
#' @param exprs optional matrix of normalized expression (observed activation scheme).
#' @param approach not used
#' @param model Model object.
#' @param bString Binary string denoting the hyper-graph.
#' @param Egenes Maximal number of visualized E-genes.
#' @param Sgene Integer denoting the S-gene. See colnames(CNOlist@signals[[1]]) to match integer with S-gene name.
#' @param parameters not used
#' @param plot Plot the heatmap. If FALSE, only corresponding information is outputed.
#' @param disc Discretize the data.
#' @param affyIds Experimental. Turn Affymetrix Ids into HGNC gene symbols.
#' @param sim not used
#' @param relFit not used
#' @param complete not used
#' @param xrot See ?heatmapOP.
#' @param Rowv See ?heatmapOP.
#' @param Colv See ?heatmapOP.
#' @param dendrogram See ?heatmapOP.
#' @param soft not used
#' @param colSideColors See ?heatmapOP.
#' @param affychip Define Affymetrix chip used to generate the data (optional).
#' @param method Scoring method can be a correlation or distance measure. See ?cor and ?dist for details.
#' @param ranks Turn data into ranks.
#' @param breaks See ?heatmapOP.
#' @param col See ?heatmapOP.
#' @param csc not used (devel)
#' @param sizeFac Size factor penelizing the hyper-graph size.
#' @param verbose Verbose output.
#' @param order Order by "rank" or "name".
#' @param colnames not used (devel)
#' @param ... additional arguments
#' @author Martin Pirkl
#' @return lattice object with matrix information
#' @export
#' @import
#' CellNOptR
#' stats
#' RColorBrewer
#' @examples
#' library(bnem)
#' library(CellNOptR)
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
#' model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnem(search = "greedy", CNOlist = CNOlist, NEMlist = NEMlist, model = model, parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
#' rownames(fc) <- 1:nrow(fc)
#' val <- validateGraph(CNOlist = CNOlist, NEMlist = NEMlist, model = model, bString = res$bString, Egenes = 10, Sgene = 4)
validateGraph <-
function(CNOlist, fc=NULL, exprs=NULL, approach = "fc", model, bString, Egenes = 25, Sgene = 1,
                          parameters = list(cutOffs = c(0,1,0), scoring = c(0.1,0.2,0.9)), plot = TRUE,
         disc = 0, affyIds = TRUE, sim = 0, relFit = FALSE, complete = FALSE, xrot = 25, Rowv = F, Colv = F, dendrogram = "none", soft = TRUE, colSideColors = NULL, affychip = "hgu133plus2", method = "s", ranks = F, breaks = NULL, col = "RdYlGn", csc = TRUE, sizeFac = 10^-10, verbose = T, order = "rank", colnames = "bio", ...) { ## order can be none, rank or names; names, rank superceed Rowv = TRUE

    NEMlist <- list()
    NEMlist$fc <- fc
    NEMlist$exprs <- exprs

  myCN2bioCN <- function(x, stimuli, inhibitors) {
    y <- gsub("_vs_", ") vs (", x)
    for (i in inhibitors) {
      y <- gsub(i, paste(i, "\\-", sep = ""), y)
    }
    for (i in stimuli) {
      y <- gsub(i, paste(i, "\\+", sep = ""), y)
    }
    y <- gsub("Ctrl", "control", paste("(", gsub("_", ",", y), ")", sep = ""))
    return(y)
  }

  gene2protein <- function(genes, strict = FALSE) {
    if (strict) {
      gene2prot <- cbind(
                     c("^APC$", "^ATF2$", "^BIRC2$", "^BIRC3$", "^CASP4$", "^CASP8$", "^CFLAR$", "^CHUK$", "^CTNNB1$", "^DKK1$", "^DKK4$", "^FLASH$", "^IKBKB$", "^IKBKG$", "^JUN$", "^MAP2K1$", "^MAP3K14$", "^MAP3K7$", "^MAPK8$", "^PIK3CA$", "^RBCK1$", "^RELA$", "^RIPK1$", "^RIPK3$", "^RNF31$", "^SHARPIN$", "^TAB2$", "^TCF4$", "^TCF7L2$", "^TNFRSF10A$", "^TNFRSF10B$", "^TNFRSF1A$", "^TNIK$", "^TRAF2$", "^USP2$", "^WLS$", "^WNT11$", "^WNT5A$", "^TNFa$", "^TRAIL"),
                     c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                     )
    } else {
      gene2prot <- cbind(
                     c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                     c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                     )
    }
    for (i in 1:nrow(gene2prot)) {
      genes <- gsub(gene2prot[i, 1], gene2prot[i, 2], genes)
    }
    return(genes)
  }

  protein2gene <- function(proteins, strict = FALSE) {
    if (strict) {
      gene2prot <- cbind(
                     c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                     c("^Apc$", "^Atf2$", "^cIap1$", "^cIap2$", "^Casp4$", "^Casp8$", "^c-Flip$", "^Ikka$", "^Beta-Cat.$", "^Dkk1$", "^Dkk4$", "^Casp8ap2$", "^Ikkb$", "^Nemo$", "^cJun$", "^Mekk$", "^Nik$", "^Tak1$", "^Jnk$", "^Pi3k$", "^Hoil1$", "^RelA$", "^Rip1$", "^Rip3$", "^Hoip$", "^Sharpin$", "^Tab2$", "^fake$", "^Tcf4$", "^Dr4$", "^Dr5$", "^Tnfr1$", "^Tnik$", "^Traf2$", "^Usp2$", "^Evi$", "^Wnt11$", "^Wnt5A$", "^Tnfa$", "^Trail$")
                     )
    } else {
      gene2prot <- cbind(
                     c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                     c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                     )
    }
    for (i in 1:nrow(gene2prot)) {
      proteins <- gsub(gene2prot[i, 2], gene2prot[i, 1], proteins)
    }
    return(proteins)
  }

  colSideColorsSave <- NULL
  bad.data <- FALSE
  errorMat <- function() {
    error.mat <- matrix(0, nrow = 50, ncol = (7*3+8))
    error.mat[1:10, c(2:4, 6:8, 10:11, 14:15, 18:20, 22:24, 26:28)] <- 3.5
    error.mat[11:20, c(2,4, 6,8, 10,12, 14,16, 18,20, 23, 26,28)] <- 3.5
    error.mat[21:30, c(2:4, 6:8, 10,12, 14,16, 18:20, 23, 26:28)] <- 3.5
    error.mat[31:40, c(2,4, 6,8, 10,12, 14,16, 18,20, 23, 26,28)] <- 3.5
    error.mat[41:50, c(2:4, 6,8, 10:11, 14:15, 18,20, 23, 26,28)] <- 3.5
    error.mat <- cbind(error.mat[, 1:13], matrix(0, nrow = 50, ncol = 5), error.mat[, 14:29])
    error.mat <- rbind(matrix(0, nrow = 10, ncol = ncol(error.mat)), error.mat, matrix(0, nrow = 10, ncol = ncol(error.mat)), error.mat, matrix(0, nrow = 10, ncol = ncol(error.mat)), error.mat, matrix(0, nrow = 10, ncol = ncol(error.mat)), error.mat, matrix(0, nrow = 10, ncol = ncol(error.mat)))
    rownames(error.mat) <- 1:250
    colnames(error.mat) <- 1:((7*3+8)+5)
    error.mat <- error.mat
    error.mat <- error.mat
    error.mat <- error.mat[rep(1:nrow(error.mat), each = 1), rep(1:ncol(error.mat), each = 5)]
    error.mat <- smoothMatrix(error.mat, 5, torus = F)
    return(error.mat)
  }
  CNOlist <- checkCNOlist(CNOlist)
  NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist, parameters = parameters, approach = approach, method = method)
  tmp <- computeScoreNemT1(CNOlist, model = model, bString, NAFac = 1, approach = approach, NEMlist = NEMlist, tellme = 1, parameters = parameters, sim = sim, relFit = relFit, method = method, sizeFac = sizeFac, verbose = F)
  EtoS <- tmp$EtoS
  if (verbose) {
    print(paste(Sgene, ".", colnames(CNOlist@signals[[1]])[Sgene], ": ", sum(EtoS[, 2] == Sgene), sep = ""))
    print(paste("Activated: ", sum(EtoS[, 2] == Sgene & EtoS[, 3] == 1), sep = ""))
    print(paste("Inhibited: ", sum(EtoS[, 2] == Sgene & EtoS[, 3] == -1), sep = ""))
    print("Summary Score:")
    print(summary(EtoS[which(EtoS[, 2] == Sgene), 4]))
    dups <- sum(duplicated(rownames(EtoS)) == TRUE)
    if (dups > 0) {
      used <- length(unique(rownames(EtoS)))
    } else {
      used <- nrow(EtoS)
    }
    print(paste("Unique genes used: ", (used), " (", round((used/nrow(NEMlist$fc))*100, 2), " %)", sep = ""))
    print(paste("Duplicated genes: ", dups, sep = ""))
    print("Overall fit:")
    print(summary(EtoS[, 4]))
  }
  indexList <- NULL
  if (is.null(indexList) == TRUE) {
    indexList = indexFinder(CNOlist, model)
  }
  modelCut = cutModel(model, bString)
  if (sim == 0) {
    simResults <- simulateStatesRecursive(CNOlist = CNOlist, model = modelCut, bString = (numeric(length(modelCut$reacID)) + 1))
  }
  if (sim == 1) {
    simResults <- simulateStatesRecursiveAdd(CNOlist = CNOlist, model = modelCut, bString = (numeric(length(modelCut$reacID)) + 1))
  }
  rownames(simResults) <- 1:nrow(simResults)
  simResults <- simResults[, which(colnames(simResults) %in% colnames(CNOlist@signals[[1]]))]
  SCompMat <- computeFc(CNOlist, t(simResults))
  SCompMat <- SCompMat[, colnames(NEMlist$fc)]

  if (parameters$cutOffs[3] == -1) {
    method <- checkMethod(method)
    S.mat <- SCompMat
    ## do median polish over gene clusters
    data.med <- NEMlist$fc[1:ncol(CNOlist@signals[[2]]), ]*0
    Epos <- EtoS[order(rownames(EtoS)), 1:2]
    for (i in 1:ncol(CNOlist@signals[[1]])) {
      tmp <- medpolish(rbind(NEMlist$fc[Epos[which(Epos[, 2] == i), 1], ], NEMlist$fc[Epos[which(Epos[, 2] == (i+ncol(CNOlist@signals[[2]]))), 1], ]), trace.iter=F)
      data.med[i, ] <- tmp$col
    }
    E.mat <- data.med
    E.mat[is.na(E.mat)] <- 0
    tmp <- which(apply(E.mat, 1, sum) != 0)
    E.mat <- E.mat[which(apply(E.mat, 1, sum) != 0), ]
    rownames(E.mat) <- rownames(S.mat)[tmp]
    NEMlist$fc <- E.mat
    if (parameters$scoring[1] == 0) {
      S.mat[which(S.mat == 0)] <- NA
    }
    if ("pearson" %in% method) {
      cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p", use = "pairwise.complete.obs"))
    }
    if ("spearman" %in% method) {
      S.mat.ranks <- apply(S.mat, 1, rank)
      cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p", use = "pairwise.complete.obs"))
    }
    cosine.sim[is.na(cosine.sim)] <- 0
    MSEAfc <- cosine.sim
    MSEIfc <- -cosine.sim
    R <- cbind(MSEAfc, MSEIfc)
    R[is.na(R)] <- max(R[!is.na(R)])
    MSEE <- rowMins(R)
  }

  ## matrix visualisation for egenes fitted:
  if ("fc" %in% approach) {
    check.data <- NEMlist$fc
    check.model <- SCompMat
    rownames(check.model) <- colnames(simResults)
  }
  if ("abs" %in% approach) {
    check.data <- NEMlist$exprs # here the $norm is needed
    check.model <- simResults
    colnames(check.model) <- colnames(CNOlist@signals[[2]])
  }

  Egenes <- Egenes

  Egenes <- min(Egenes, sum(EtoS[, 2] == Sgene))

  if (Egenes == 0) {
    mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
    if (plot) {
      print(heatmapOP(errorMat(), main = mainlab, sub = "", Colv = F, Rowv = F, col = "RdBu", coln = 12, breaks = seq(0,8, 0.1)), ...)
    }
    genesInfo <- as.matrix(0)
    rownames(genesInfo) <- "dummy"
    return(list(genesInfo = genesInfo, data = genesInfo))
  }

  if (complete) {
    egenefit <- matrix(0, nrow = (sum(EtoS[, 2] == Sgene)+1), ncol = ncol(check.data))
  } else {
    egenefit <- matrix(0, nrow = (Egenes+1), ncol = ncol(check.data))
  }

  egenefit[1,] <- check.model[Sgene, ]
  rownames(egenefit) <- 1:nrow(egenefit)
  rownames(egenefit)[1] <- rownames(check.model)[Sgene]
  colnames(egenefit) <- 1:ncol(egenefit)
  colnames(egenefit) <- colnames(check.data)
  #colnames(egenefit) <- gsub(".*_vs_", "", gsub("Ctrl_vs_", "", colnames(check.data)))
  if (complete) {
    activatedEgenes <- numeric(sum(EtoS[, 2] == Sgene)+1)
  } else {
    activatedEgenes <- numeric(Egenes+1)
  }

  count <- 0
  for (i in 1:nrow(EtoS)) {
    if (EtoS[i, 2] == Sgene) {
      egenefit[count+2,] <- check.data[which(rownames(check.data) %in% rownames(EtoS)[i]), ]
      rownames(egenefit)[count+2] <- rownames(EtoS)[i]
      if (EtoS[i, 3] == 1) {
        activatedEgenes[count+2] <- 1
      }
      count <- count + 1
    }
    if (count >= Egenes & !complete) { break() }
  }

  Egenes <- count

  if (affyIds == FALSE) {
    require(paste(affychip, ".db", sep = ""), character.only = T) # is this ok?
    temp <- as.vector(unlist(mget(unique(rownames(egenefit)[-1]), get(paste(affychip, "SYMBOL", sep = "")))))
    temp2 <- rownames(egenefit)[-1]
    if (sum(is.na(temp) == T) > 0) {
      temp[is.na(temp)] <- temp2[is.na(temp)]
    }
    rownames(egenefit)[-1] <- temp
  }

  rownames(egenefit)[is.na(rownames(egenefit))] <- "NA"
  count <- 0
  if (min(egenefit) != max(egenefit)) {
    if (disc == 1) {
      egenefit[1, ] <- disc(egenefit[1, ], 0.5)
      for (i in 2:nrow(egenefit)) {
        if (parameters$cutOffs[2] == 0) {
          save.model <- egenefit[1, ]
          egenefit[1, ] <- save.model
          real.breaks <- seq(-max(abs(min(egenefit)),abs(max(egenefit))),max(abs(min(egenefit)),abs(max(egenefit))),0.1)
          if (length(real.breaks) > 101) {
            real.breaks <- real.breaks[c((floor(length(real.breaks)/2)-50):floor(length(real.breaks)/2), (floor(length(real.breaks)/2)+1):(floor(length(real.breaks)/2)+50))]
          }
        } else {
          if (activatedEgenes[i] == 1) {
            onePosI <- c(which(egenefit[1, ] == 1 & egenefit[i, ] > parameters$cutOffs[2]),
                         which(egenefit[1, ] == 0 & egenefit[i, ] > parameters$cutOffs[2]),
                         which(egenefit[1, ] == -1 & egenefit[i, ] > parameters$cutOffs[2]))
            onePosII <- which(egenefit[1, ] == 1 & egenefit[i, ] > parameters$cutOffs[1] & egenefit[i, ] <= parameters$cutOffs[2])
            oneNegI <- c(which(egenefit[1, ] == -1 & egenefit[i, ] < -parameters$cutOffs[2]),
                         which(egenefit[1, ] == 0 & egenefit[i, ] < -parameters$cutOffs[2]),
                         which(egenefit[1, ] == 1 & egenefit[i, ] < -parameters$cutOffs[2]))
            oneNegII <- which(egenefit[1, ] == -1 & egenefit[i, ] < -parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[2])
            zeros <- c(which(egenefit[1, ] == 0 & abs(egenefit[i, ]) <= parameters$cutOffs[2]),
                       which(egenefit[1, ] == 1 & egenefit[i, ] <= parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[2]),
                       which(egenefit[1, ] == -1 & egenefit[i, ] >= -parameters$cutOffs[1] & egenefit[i, ] <= parameters$cutOffs[2]))
            zerosI <- c(which(egenefit[1, ] == 0 & abs(egenefit[i, ]) <= parameters$cutOffs[1]),
                        which(egenefit[1, ] == 1 & egenefit[i, ] <= parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[1]),
                        which(egenefit[1, ] == -1 & egenefit[i, ] <= parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[1]))
            zerosII <- c(which(egenefit[1, ] == 0 & egenefit[i, ] <= parameters$cutOffs[2] & egenefit[i, ] > parameters$cutOffs[1]),
                         which(egenefit[1, ] == -1 & egenefit[i, ] > parameters$cutOffs[1] & egenefit[i, ] <= parameters$cutOffs[2]))
            zerosIII <- c(which(egenefit[1, ] == 0 & egenefit[i, ] >= -parameters$cutOffs[2] & egenefit[i, ] < -parameters$cutOffs[1]),
                          which(egenefit[1, ] == 1 & egenefit[i, ] < -parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[2]))
            if (soft) {
              egenefit[i, onePosI] <- 3
              egenefit[i, oneNegI] <- -3
              egenefit[i, onePosII] <- 2
              egenefit[i, oneNegII] <- -2
              egenefit[i, zerosI] <- 0
              egenefit[i, zerosII] <- 0.5
              egenefit[i, zerosIII] <- -0.5
            } else {
              egenefit[i, onePosI] <- 2
              egenefit[i, oneNegI] <- -2
              egenefit[i, onePosII] <- 2
              egenefit[i, oneNegII] <- -2
              egenefit[i, zerosI] <- 0
              egenefit[i, zerosII] <- 0
              egenefit[i, zerosIII] <- 0
            }
          } else {
            onePosI <- c(which(egenefit[1, ] == 1 & egenefit[i, ] > parameters$cutOffs[2]),
                         which(egenefit[1, ] == 0 & egenefit[i, ] > parameters$cutOffs[2]),
                         which(egenefit[1, ] == -1 & egenefit[i, ] > parameters$cutOffs[2]))
            onePosII <- which(egenefit[1, ] == -1 & egenefit[i, ] > parameters$cutOffs[1] & egenefit[i, ] <= parameters$cutOffs[2])
            oneNegI <- c(which(egenefit[1, ] == -1 & egenefit[i, ] < -parameters$cutOffs[2]),
                         which(egenefit[1, ] == 0 & egenefit[i, ] < -parameters$cutOffs[2]),
                         which(egenefit[1, ] == 1 & egenefit[i, ] < -parameters$cutOffs[2]))
            oneNegII <- which(egenefit[1, ] == 1 & egenefit[i, ] < -parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[2])
            zerosI <- c(which(egenefit[1, ] == 0 & abs(egenefit[i, ]) <= parameters$cutOffs[1]),
                        which(egenefit[1, ] == 1 & egenefit[i, ] <= parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[1]),
                        which(egenefit[1, ] == -1 & egenefit[i, ] <= parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[1]))
            zerosII <- c(which(egenefit[1, ] == 0 & egenefit[i, ] <= parameters$cutOffs[2] & egenefit[i, ] > parameters$cutOffs[1]),
                         which(egenefit[1, ] == 1 & egenefit[i, ] > parameters$cutOffs[1] & egenefit[i, ] <= parameters$cutOffs[2]))
            zerosIII <- c(which(egenefit[1, ] == 0 & egenefit[i, ] >= -parameters$cutOffs[2] & egenefit[i, ] < -parameters$cutOffs[1]),
                          which(egenefit[1, ] == -1 & egenefit[i, ] < -parameters$cutOffs[1] & egenefit[i, ] >= -parameters$cutOffs[2]))
            if (soft) {
              egenefit[i, onePosI] <- 3
              egenefit[i, oneNegI] <- -3
              egenefit[i, onePosII] <- 2
              egenefit[i, oneNegII] <- -2
              egenefit[i, zerosI] <- 0
              egenefit[i, zerosII] <- 0.5
              egenefit[i, zerosIII] <- -0.5
            } else {
              egenefit[i, onePosI] <- 2
              egenefit[i, oneNegI] <- -2
              egenefit[i, onePosII] <- 2
              egenefit[i, oneNegII] <- -2
              egenefit[i, zerosI] <- 0
              egenefit[i, zerosII] <- 0
              egenefit[i, zerosIII] <- 0
            }
          }
          real.breaks <- c(-4.5,-3.5,-2.5,-1.5,-0.75,-0.25,0.25,0.75,1.5,2.5,3.5,4.5)
        }
      }
      egenefit[1, ] <- egenefit[1, ]*4
      if (!is.null(colSideColors)) {
        colSideColorsSave <- colSideColors
      }
      if (min(egenefit) != max(egenefit)) {
        if (csc) {
          colSideColors <- rep("black", ncol(egenefit))
          for (i in 1:length(colnames(egenefit))) {
            names <- unlist(strsplit(colnames(egenefit)[i], "_"))
            if (length(names) > 1) {
              if (names[1] %in% colnames(CNOlist@inhibitors) & names[length(names)] %in% colnames(CNOlist@stimuli)) {
                colSideColors[i] <- "brown"
              }
              if (names[1] %in% colnames(CNOlist@stimuli) & names[length(names)] %in% colnames(CNOlist@inhibitors)) {
                colSideColors[i] <- "orange"
              }
              if (names[1] %in% "Ctrl" & names[length(names)] %in% colnames(CNOlist@stimuli)) {
                colSideColors[i] <- "yellow"
              }
              if (names[1] %in% "Ctrl" & names[length(names)] %in% colnames(CNOlist@inhibitors)) {
                colSideColors[i] <- "blue"
              }
            } else {
              if (names %in% colnames(CNOlist@inhibitors)) {
                colSideColors[i] <- "blue"
              } else {
                colSideColors[i] <- "yellow"
              }
            }
          }
          if (!is.null(colSideColorsSave)) {
            colSideColors <- rbind(colSideColorsSave, colSideColors)
          }
          order2 <- order(colSideColors)
        }
        order1 <- order(colnames(egenefit))
        egenefit <- egenefit[nrow(egenefit):1, ]
        if (nrow(egenefit) == 2) {
          names.backup <- rownames(egenefit)[1]
          egenefit <- rbind(egenefit[1:(nrow(egenefit)-1), ], 0, egenefit[nrow(egenefit), ])
          rownames(egenefit) <- 1:nrow(egenefit)
          rownames(egenefit)[1] <- names.backup
        } else {
          egenefit <- rbind(egenefit[1:(nrow(egenefit)-1), ], 0, egenefit[nrow(egenefit), ])
        }
        rownames(egenefit)[nrow(egenefit)] <- rownames(check.model)[Sgene]
        rownames(egenefit)[(nrow(egenefit)-1)] <- "---"
        if (parameters$scoring[1] == 0) {
          ##col.sums <- apply(abs(egenefit), 2, median)
          col.sums <- colMedians(abs(egenefit))
          sig.mismatch <- which(col.sums >= 2)
          sig.mismatch <- sig.mismatch[-which(egenefit[nrow(egenefit), ] != 0)]
          get.cols <- unique(c(sig.mismatch, which(egenefit[nrow(egenefit), ] != 0)))
          egenefit <- egenefit[, get.cols]
        } else {
          get.cols <- 1:ncol(egenefit)
        }
        mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
        if (csc) {
          colSideColors <- colSideColors[get.cols]
          if ("yellow" %in% colSideColors) {
            mainlab <- paste(mainlab, " Stimulationeffects (yellow)", sep = "")
          }
          if ("blue" %in% colSideColors) {
            mainlab <- paste(mainlab, " Knockdowneffects (blue)\n", sep = "")
          }
          if ("orange" %in% colSideColors) {
            mainlab <- paste(mainlab, " Silencingeffects Type I (orange)", sep = "")
          }
          if ("brown" %in% colSideColors) {
            mainlab <- paste(mainlab, " Silencingeffects Type II (brown)", sep = "")
          }
        } else {
          colSideColors <- NULL
        }
        if (order %in% "rank") {
          Rowv <- FALSE
          egenefit_genes <- egenefit[1:(nrow(egenefit)-2), ]
          namereset <- FALSE
          if (!is.matrix(egenefit_genes)) {
            egenefit_genes <- t(as.matrix(egenefit_genes))
            namereset <- TRUE
          } else {
            geneorder <- rownames(EtoS)[which(EtoS[, 2] == Sgene)[1:Egenes]]
            egenefit_genes <- egenefit_genes[order(match(rownames(egenefit_genes), geneorder), decreasing = T), ]
          }
          egenefit <- rbind(egenefit_genes, egenefit[(nrow(egenefit)-1):nrow(egenefit), ])
          if (namereset) {
            rownames(egenefit)[1] <- names(which(EtoS[, 2] == Sgene))
          }
        }
        if (order %in% "names") {
          Rowv <- FALSE
          tmp <- egenefit[1:(nrow(egenefit)-2), ]
          tmp <- tmp[order(rownames(tmp)), ]
          rownames(egenefit)[1:(nrow(egenefit)-2)] <- rownames(tmp)
          egenefit[1:(nrow(egenefit)-2), ] <- tmp
        }
        if (!is.null(dim(egenefit)) & plot) {
          if (Rowv & nrow(egenefit) > 3) {
            Rowv <- F
            d <- dist(egenefit[-c(nrow(egenefit)-1, nrow(egenefit)), ])
            hc <- hclust(d)
            egenefit <- egenefit[c(hc$order, nrow(egenefit)-1, nrow(egenefit)), ]
          } else {
            Rowv = F
          }
          if (Colv) {
            clusterdata <- egenefit
            clusterdata[nrow(egenefit), ] <- 0
            clusterdata <- NULL
          }
          clusterdata <- NULL
          if ("bio" %in% colnames) {
            colnames(egenefit) <- gene2protein(myCN2bioCN(colnames(egenefit), colnames(CNOlist@stimuli), colnames(CNOlist@inhibitors)))
          }
          print(heatmapOP(egenefit, main = mainlab, xrot = xrot, breaks = real.breaks, coln = 11, Colv = Colv, Rowv = Rowv, colSideColors = colSideColors, dendrogram = dendrogram, col = col, clusterx = clusterdata, ...))
        } else {
          print("one effect is not a matrix")
        }
      } else {
        print("min equals max in data matrix")
      }
    } else {
      if (csc) {
        colSideColors <- rep("black", ncol(egenefit))
        for (i in 1:length(colnames(egenefit))) {
          names <- unlist(strsplit(colnames(egenefit)[i], "_"))
          if (length(names) > 1) {
            if (names[1] %in% colnames(CNOlist@inhibitors) & names[length(names)] %in% colnames(CNOlist@stimuli)) {
              colSideColors[i] <- "brown"
            }
            if (names[1] %in% colnames(CNOlist@stimuli) & names[length(names)] %in% colnames(CNOlist@inhibitors)) {
              colSideColors[i] <- "orange"
            }
            if (names[1] %in% "Ctrl" & names[length(names)] %in% colnames(CNOlist@stimuli)) {
              colSideColors[i] <- "yellow"
            }
            if (names[1] %in% "Ctrl" & names[length(names)] %in% colnames(CNOlist@inhibitors)) {
              colSideColors[i] <- "blue"
            }
          } else {
            if (names %in% colnames(CNOlist@inhibitors)) {
              colSideColors[i] <- "blue"
            } else {
              colSideColors[i] <- "yellow"
            }
          }
        }
        if (!is.null(colSideColorsSave)) {
          colSideColors <- rbind(colSideColorsSave, colSideColors)
        }
        order2 <- order(colSideColors)
      }
      order1 <- order(colnames(egenefit))
      egenefit <- egenefit[nrow(egenefit):1, ]
      if (nrow(egenefit) == 2) {
        names.backup <- rownames(egenefit)[1]
        egenefit <- rbind(egenefit[1:(nrow(egenefit)-1), ], 0, egenefit[nrow(egenefit), ])
        rownames(egenefit) <- 1:nrow(egenefit)
        rownames(egenefit)[1] <- names.backup
      } else {
        egenefit <- rbind(egenefit[1:(nrow(egenefit)-1), ], 0, egenefit[nrow(egenefit), ])
      }
      rownames(egenefit)[nrow(egenefit)] <- rownames(check.model)[Sgene]
      rownames(egenefit)[(nrow(egenefit)-1)] <- "---"
      if (parameters$scoring[1] == 0) {
        ##col.sums <- apply(abs(egenefit), 2, median)
        col.sums <- colMedians(abs(egenefit))
        sig.mismatch <- which(col.sums >= 2)
        sig.mismatch <- sig.mismatch[-which(egenefit[nrow(egenefit), ] != 0)]
        get.cols <- unique(c(sig.mismatch, which(egenefit[nrow(egenefit), ] != 0)))
        egenefit <- egenefit[, get.cols]
      } else {
        get.cols <- 1:ncol(egenefit)
      }
      mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
      if (csc) {
        colSideColors <- colSideColors[get.cols]
        if ("yellow" %in% colSideColors) {
          mainlab <- paste(mainlab, " Stimulationeffects (yellow)", sep = "")
        }
        if ("blue" %in% colSideColors) {
          mainlab <- paste(mainlab, " Knockdowneffects (blue)\n", sep = "")
        }
        if ("orange" %in% colSideColors) {
          mainlab <- paste(mainlab, " Silencingeffects Type I (orange)", sep = "")
        }
        if ("brown" %in% colSideColors) {
          mainlab <- paste(mainlab, " Silencingeffects Type II (brown)", sep = "")
        }
      } else {
        colSideColors <- NULL
      }
      mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
      egenefit[nrow(egenefit), ] <- egenefit[nrow(egenefit), ]*max(abs(egenefit))
      if (order %in% "rank") {
        Rowv <- FALSE
        egenefit_genes <- egenefit[1:(nrow(egenefit)-2), ]
        namereset <- FALSE
        if (!is.matrix(egenefit_genes)) {
          egenefit_genes <- t(as.matrix(egenefit_genes))
          namereset <- TRUE
        } else {
          geneorder <- rownames(EtoS)[which(EtoS[, 2] == Sgene)[1:Egenes]]
          egenefit_genes <- egenefit_genes[order(match(rownames(egenefit_genes), geneorder), decreasing = T), ]
        }
        egenefit <- rbind(egenefit_genes, egenefit[(nrow(egenefit)-1):nrow(egenefit), ])
        if (namereset) {
          rownames(egenefit)[1] <- names(which(EtoS[, 2] == Sgene))
        }
      }
      if (order %in% "names") {
        Rowv <- FALSE
        tmp <- egenefit[1:(nrow(egenefit)-2), ]
        tmp <- tmp[sort(rownames(tmp)), ]
        rownames(egenefit)[1:(nrow(egenefit)-2)] <- rownames(tmp)
        egenefit[1:(nrow(egenefit)-2), ] <- tmp
      }
      if (Rowv & nrow(egenefit) > 3) {
        Rowv <- F
        d <- dist(egenefit[-c(nrow(egenefit)-1, nrow(egenefit)), ])
        hc <- hclust(d)
        egenefit <- egenefit[c(hc$order, nrow(egenefit)-1, nrow(egenefit)), ]
      } else {
        Rowv = F
      }
      if (Colv) {
        clusterdata <- egenefit
        clusterdata[nrow(egenefit), ] <- 0
        clusterdata <- NULL
      }
      clusterdata <- NULL
      low <- sum(egenefit[nrow(egenefit), ] == min(egenefit[nrow(egenefit), ]))
      high <- sum(egenefit[nrow(egenefit), ] == max(egenefit[nrow(egenefit), ]))
      egenefit2 <- egenefit
      egenefit2[which(rownames(egenefit2) %in% rownames(EtoS)[which(EtoS[, 3] == -1)]), ] <- egenefit2[which(rownames(egenefit2) %in% rownames(EtoS)[which(EtoS[, 3] == -1)]), ]*(-1)
      egenefit2 <- t(apply(egenefit2, 1, rank))
      high2 <- which(egenefit2 > ncol(egenefit2)-high)
      low2 <- which(egenefit2 < low)
      mid <- which(egenefit2 >= low & egenefit2 <= ncol(egenefit2)-high)
      egenefit2[high2] <- 1
      egenefit2[low2] <- -1
      egenefit2[mid] <- 0
      if (csc) {
        colSideColors <- colSideColors[order(egenefit2[nrow(egenefit2), ])]
      }
      egenefit <- egenefit[, order(egenefit2[nrow(egenefit2), ])]
      egenefit2 <- egenefit2[, order(egenefit2[nrow(egenefit2), ])]
      show.ranks <- FALSE#T
      if (show.ranks) {
        egenefit <- t(apply(egenefit, 1, rank))
        breaks <- c(0, sort(unique(egenefit[nrow(egenefit), ])), ncol(egenefit)+1)
        print(breaks)
      }
      if (ranks) {
        ##egenefit2 <- rbind(egenefit2[1:(nrow(egenefit2)-2), ], egenefit2[(nrow(egenefit2)-1):nrow(egenefit2), ])
        if (is.null(breaks)) {
          breaks <- c(-2,-0.5,0.5,2)
        }
        if ("bio" %in% colnames) {
          colnames(egenefit2) <- gene2protein(myCN2bioCN(colnames(egenefit2), colnames(CNOlist@stimuli), colnames(CNOlist@inhibitors)))
        }
        print(heatmapOP(egenefit2, main = mainlab, xrot = xrot, breaks = breaks, coln = 11, Colv = Colv, Rowv = Rowv, colSideColors = colSideColors, dendrogram = dendrogram, colSideColorsPos = "top", col = col, clusterx = clusterdata, ...))
        ##print(heatmapOP(egenefit, main = mainlab, xrot = xrot, breaks = c(0,low,ncol(egenefit)-high,ncol(egenefit)), coln = 11, Colv = Colv, Rowv = Rowv, colSideColors = colSideColors, dendrogram = dendrogram, ...))
        ##print(heatmapOP(egenefit, main = mainlab, xrot = xrot, breaks = seq(-2,2,0.1), coln = 11, Colv = Colv, Rowv = Rowv, colSideColors = colSideColors, dendrogram = dendrogram, ...))
      } else {
        ##egenefit <- rbind(egenefit[nrow(egenefit):(nrow(egenefit)-1), ], egenefit[1:(nrow(egenefit)-2), ])
        if (is.null(breaks)) {
          breaks <- seq(-1,1,0.1)
        }
        if ("bio" %in% colnames) {
          colnames(egenefit) <- gene2protein(myCN2bioCN(colnames(egenefit), colnames(CNOlist@stimuli), colnames(CNOlist@inhibitors)))
        }
        print(heatmapOP(egenefit, main = mainlab, xrot = xrot, breaks = breaks, coln = 11, Colv = Colv, Rowv = Rowv, colSideColors = colSideColors, dendrogram = dendrogram, colSideColorsPos = "top", col = col, clusterx = clusterdata, ...))
      }
    }
  } else {
    mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
    if (plot) {
      print(heatmapOP(errorMat(), main = mainlab, sub = "", Colv = F, Rowv = F, col = "RdBu", coln = 12, breaks = seq(0,8,0.1)), ...)
    }
    print("min equals max in data matrix")
    bad.data <- TRUE
  }
  if (sum(EtoS[, 2] == Sgene) == 0) {
    if (!bad.data) {
      mainlab <- paste("Regulated by ", rownames(check.model)[Sgene], "\n", sep = "")
      if (plot) {
        print(heatmapOP(errorMat(), main = mainlab, sub = "", Colv = F, Rowv = F, col = "RdBu", coln = 12, breaks = seq(0,8,0.1)), ...)
      }
      return(NULL)
    }
  } else {
    if (!is.na(rownames(as.matrix(EtoS[which(EtoS[, 2] == Sgene), ]))[1])) {
      if (sum(EtoS[, 2] == Sgene) > 1) {
        genesInfo <- EtoS[which(EtoS[, 2] == Sgene), ]
      } else {
        names.backup <- rownames(EtoS)[which(EtoS[, 2] == Sgene)]
        genesInfo <- t(as.matrix(EtoS[which(EtoS[, 2] == Sgene), ]))
        rownames(genesInfo) <- names.backup
      }
      if (affyIds == FALSE) {
        require(paste(affychip, ".db", sep = ""), character.only = T)
        temp <- as.vector(unlist(mget(rownames(genesInfo), get(paste(affychip, "SYMBOL", sep = "")))))
        if (sum(is.na(temp) == T) > 0) {
          temp[is.na(temp)] <- rownames(genesInfo)[is.na(temp)]
        }
        rownames(genesInfo) <- temp
      }
      return(list(genesInfo = genesInfo, data = egenefit))
    }
  }
}
