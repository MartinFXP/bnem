crossTalk <- function(CNOlist,
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
                      Sgenes = NULL
                      ) {
  stimuli <- colnames(CNOlist@stimuli)
  inhibitorsFull <- colnames(CNOlist@inhibitors)
  if (is.null(Sgenes)) {
    Sgenes <- inhibitorsFull
  }
  allCross <- list()
  for (i in inhibitorsFull) {
    if (!(i %in% Sgenes)) { next() }
    CNOlistTmp <- CNOlist
    inhibitors <- i
    CNOlistTmp@cues <- CNOlist@cues[, c(stimuli, i)]
    CNOlistTmp@signals[[1]] <- as.matrix(CNOlist@signals[[1]][, i])
    colnames(CNOlistTmp@signals[[1]]) <- i
    CNOlistTmp@signals[[2]] <- as.matrix(CNOlist@signals[[2]][, i])
    colnames(CNOlistTmp@signals[[2]]) <- i
    CNOlistTmp@signals[[1]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[1]])
    CNOlistTmp@signals[[2]] <- cbind(CNOlistTmp@stimuli, CNOlistTmp@signals[[2]])
    CNOlistTmp@stimuli <- CNOlist@stimuli[, stimuli]
    CNOlistTmp@inhibitors <- as.matrix(CNOlist@inhibitors[, i])
    colnames(CNOlistTmp@inhibitors) <- i
    
    sifMatrix <- numeric()
    signals <- c(stimuli, inhibitors)
    for (j in stimuli) {
      if (i %in% j) { next() }
      sifMatrix <- rbind(sifMatrix, c(j, 1, i))
      sifMatrix <- rbind(sifMatrix, c(j, -1, i))
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
    cnoStimuli <- apply(CNOlist@stimuli[, stimuli], 1, sum)
    cnoInhibitors <- CNOlist@inhibitors[, i]
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
