triples <- function(CNOlist,
                    initBstring = NULL, # initBstring = TRUE
                    sizeFac = 1e-04, 
                    NAFac = 1,
                    popSize = 50,
                    pMutation = 0.5,
                    maxTime = 60,
                    maxGens = 500, 
                    stallGenMax = 100,
                    relTol = 0.1, 
                    verbose = TRUE,
                    priorBitString = NULL,
                    selPress = c(1.2,0.0001), # 1.2
                    approach = "fc",
                    NEMlist,
                    fit = "linear",
                    targetBstring = "none",
                    elitism = NULL,
                    inversion = NULL,
                    graph = FALSE,
                    parameters = list(cutOffs = c(0.7,0.7,0.7), scoring = c(2,1), fitBest = 0),
                    parallel = NULL, # parallel = N with N number of cores to use or a list with cores in the first and machines in the second entry like list(cores=c(2,4,8), machines=c("bionform1", "bioinform2", "bioinform3"))
                    parallel2 = 1,
                    selection = c("tournament", 2), # can be c("tournament", k) or "sus"
                    Sgenes = NULL,
                    active = NULL
                    ) {
  stimuli <- colnames(CNOlist@stimuli)
  inhibitorsFull <- colnames(CNOlist@inhibitors)
  if (is.null(Sgenes)) {
    triples <- expand.grid(inhibitorsFull, inhibitorsFull, inhibitorsFull)
  } else {
    triples <- expand.grid(Sgenes[[1]], Sgenes[[2]], Sgenes[[3]])
  }
  allCross <- list()
  triples <- apply(triples, c(1,2), as.character)
  doubles <- NULL
  for (i in 1:nrow(triples)) {
    triples[i, ] <- sort(triples[i, ])
    if (length(unique(triples[i, ])) != length(triples[i, ])) {
      doubles <- c(doubles, i)
    }
  }
  if (length(doubles) > 0) {
    triples <- triples[-doubles, ]
  }
  triples.con <- apply(triples, 1, paste, sep = "_")
  triples.dups <- which(duplicated(triples))
  if (length(triples.dups) > 0) {
    triples <- triples[-triples.dups, ]
  }
  for (i in 1:nrow(triples)) {
    CNOlistTmp <- CNOlist
    inhibitors <- triples[i, ]
    CNOlistTmp@cues <- CNOlist@cues[, which(colnames(CNOlist@cues) %in% c(stimuli, inhibitors))]
    CNOlistTmp@signals[[1]] <- as.matrix(CNOlist@signals[[1]][, which(colnames(CNOlist@signals[[1]]) %in% inhibitors)])
    colnames(CNOlistTmp@signals[[1]]) <- inhibitors
    CNOlistTmp@signals[[2]] <- as.matrix(CNOlist@signals[[2]][, which(colnames(CNOlist@signals[[1]]) %in% inhibitors)])
    colnames(CNOlistTmp@signals[[2]]) <- inhibitors
    if (is.null(active)) {
      CNOlistTmp@signals[[1]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[1]])
      CNOlistTmp@signals[[2]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[2]])
    } else {
      CNOlistTmp@signals[[1]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[1]], 0)
      CNOlistTmp@signals[[2]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[2]], 0)
      colnames(CNOlistTmp@signals[[1]])[ncol(CNOlistTmp@signals[[1]])] <- "FALSE"
      colnames(CNOlistTmp@signals[[2]])[ncol(CNOlistTmp@signals[[2]])] <- "FALSE"
    }
    CNOlistTmp@stimuli <- CNOlist@stimuli[, which(colnames(CNOlist@stimuli) %in% stimuli)]
    CNOlistTmp@inhibitors <- as.matrix(CNOlist@inhibitors[, which(colnames(CNOlist@inhibitors) %in% inhibitors)])
    colnames(CNOlistTmp@inhibitors) <- inhibitors
    
    sifMatrix <- numeric()
    signals <- c(stimuli, inhibitors)
    for (j in stimuli) {
      ##if (i %in% j) { next() }
      sifMatrix <- rbind(sifMatrix, c(j, 1, inhibitors))
      ##sifMatrix <- rbind(sifMatrix, c(j, -1, inhibitors))
    }
    for (j in inhibitors) {
      for (k in inhibitors) {
        if (j %in% k) { next() }
        sifMatrix <- rbind(sifMatrix, c(k, 1, j, "", ""))
        ##sifMatrix <- rbind(sifMatrix, c(k, -1, j, "", ""))
      }
    }
    if (!is.null(active)) {
      for (j in active) {
        sifMatrix <- rbind(sifMatrix, c("FALSE", -1, j, "", ""))
      }
    }
    write.table(sifMatrix, file = "cstmp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    PKN <- readSIF("cstmp.sif")
    unlink("cstmp.sif")
    
    checkSignals(CNOlistTmp,PKN)
    indices<-indexFinder(CNOlistTmp,PKN,verbose=TRUE)
    NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
    NCNOcut<-cutNONC(PKN,NCNOindices)
    indicesNCNOcut<-indexFinder(CNOlistTmp,NCNOcut)
    NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
    indicesNCNOcutComp<-indexFinder(CNOlistTmp,NCNOcutComp)
    NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate=100)
    resECNOlist<-residualError(CNOlistTmp)
    Fields4Sim<-prep4sim(NCNOcutCompExp)
    model <- NCNOcutCompExp
    
    initBstring<-rep(0,length(NCNOcutCompExp$reacID))
    
    NEMlistTmp <- NEMlist
    ## the following reduction of the set is not needed! But should be used to decrease computation time for exhaustive AND genetic search:
    cnoStimuli <- apply(CNOlist@stimuli[, which(colnames(CNOlist@stimuli) %in% stimuli)], 1, sum)
    cnoInhibitors <- apply(CNOlist@inhibitors[, which(colnames(CNOlist@inhibitors) %in% inhibitors)], 1, sum)
    cnoInhibitors2 <- apply(CNOlist@inhibitors, 1, sum)
    cnoCues <- apply(CNOlist@cues, 1, sum)
    
    grepCtrls <- which(cnoCues == 0)
    grepStimuli <- which(cnoStimuli != 0 & cnoInhibitors2 == 0)
    grepInhibitors <- which(cnoStimuli == 0 & cnoInhibitors != 0)
    grepComb <- which(cnoStimuli != 0 & cnoInhibitors != 0)
    
    NEMlistTmp$exprs <- NEMlist$exprs[, sort(c(grepCtrls, grepStimuli, grepInhibitors, grepComb))]
    CNOlistTmp <- cnoFromData(NEMlistTmp$exprs, times = CNOlist@timepoints, stimuli = stimuli, inhibitors = inhibitors, signals = signals)
    grepFcs <- colnames(computeFc(CNOlistTmp, NEMlistTmp$exprs))
    NEMlistTmp$fc <- NEMlistTmp$fc[, which(colnames(NEMlist$fc) %in% grepFcs)]
    NEMlistTmp$B <- NEMlistTmp$B[, which(colnames(NEMlist$B) %in% grepFcs)]
    NEMlistTmp$pvals <- NEMlistTmp$pvals[, which(colnames(NEMlist$pvals) %in% grepFcs)]

    #if (parameters$cutOffs[1] != 0) {
    #  genes.max <- apply(abs(NEMlistTmp$fc), 1, max)
    #  NEMlistTmp$fc <- NEMlistTmp$fc[which(genes.max >= parameters$cutOffs[1]), ]
    #}
    
    res <- gaBinaryNemT1(parallel=parallel,
                                          CNOlist=CNOlistTmp,
                                          NEMlist=NEMlistTmp,
                                          model=model,
                                          initBstring=initBstring,
                                          popSize = popSize, 
                                          maxTime = maxTime,
                                          maxGens = maxGens,
                                          stallGenMax = stallGenMax,
                                          elitism = ceiling(popSize*0.1),
                                          inversion = ceiling(popSize*0.1),
                                          verbose=verbose,
                                          selection = selection,
                                          parameters = parameters,
                                          fit=fit,
                                          relTol=relTol,
                                          approach=approach
                                          )


    allCross[[i]] <- list(result = res, model = model, CNOlist = CNOlistTmp, NEMlist = NEMlistTmp)
    
  }
  return(allCross)
}
