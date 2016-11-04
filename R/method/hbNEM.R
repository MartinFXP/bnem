hbNEM <- function (CNOlist, PKN,
                           model,
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
                           parameters = list(cutOffs = c(0.7,0.7,0), scoring = c(2,0.5), fitBest = 0),
                           parallel = NULL, # parallel = N with N number of cores to use or a list with cores in the first and machines in the second entry like list(cores=c(2,4,8), machines=c("bionform1", "bioinform2", "bioinform3"))
                           parallel2 = 1,
                   selection = "sus" # can also be c("tournament", k)
                   ) {
  
  times <- CNOlist@timepoints
  
  if (length(initBstring) > 1) {
    if (is.matrix(initBstring)) {
    } else {
      initGates <- model$reacID[which(initBstring == 1)]
    }
  } else {
    initGates <- NULL
  }
  pGraph <- partitionGraph(model)
  stimuli <- colnames(CNOlist@stimuli)
  inhibitors <- c(colnames(CNOlist@inhibitors), model$namesSpecies[-which(model$namesSpecies %in% c(stimuli,colnames(CNOlist@inhibitors)))])
  #NEMlist$signalStates <- matrix(NA, nrow = ncol(NEMlist$exprs), ncol = length(inhibitors))
  stimuliStates <- CNOlist@stimuli
  #rownames(NEMlist$signalStates) <- colnames(NEMlist$exprs)
  #colnames(NEMlist$signalStates) <- inhibitors
  repFac <- nrow(NEMlist$signalStates)/nrow(stimuliStates)
  #NEMlist$signalStates <- cbind(stimuliStates[rep(1:nrow(stimuliStates), repFac), ], NEMlist$signalStates)
  graphs <- list()
  for (i in 1:length(pGraph$vertices)) {

    ## reduce data set to ensure consistency with subgraph
    NEMlistTmp <- NEMlist
    stimuli <- colnames(CNOlist@stimuli)
    if (sum(pGraph$vertices[[i]] %in% stimuli) > 0) {
      inhibitors <- pGraph$vertices[[i]][-which(pGraph$vertices[[i]] %in% stimuli)]
    } else {
      inhibitors <- pGraph$vertices[[i]]
    }
    ## the following reduction of the set is not needed! But should be used to decrease computation time for exhaustive AND genetic search:
    cnoStimuli <- apply(CNOlist@stimuli, 1, sum)
    cnoInhibitors <- apply(CNOlist@inhibitors, 1, sum)
    cnoCues <- apply(CNOlist@cues, 1, sum)
    grepCtrls <- which(cnoCues == 0)
    grepStimuli <- which(cnoStimuli != 0 & cnoInhibitors == 0)
    grepInhibitors <- which(cnoStimuli == 0 & cnoInhibitors != 0)
    if (sum(pGraph$vertices[[i]] %in% stimuli) > 0) {
      vertices <- pGraph$vertices[[i]][-which(pGraph$vertices[[i]] %in% stimuli)]
    } else {
      vertices <- pGraph$vertices[[i]]
    }
    NEMlistTmp$exprs <- NEMlistTmp$exprs[, sort(c(grepCtrls, grepStimuli, grep(paste(vertices, collapse = "|"), colnames(NEMlist$exprs))))]
    CNOlistTmp <- cnoFromData(NEMlistTmp$exprs, times = times, stimuli = stimuli, inhibitors = inhibitors, signals = signals)
    CNOlistTmp@cues <- CNOlistTmp@cues[, unique(colnames(CNOlistTmp@cues))]
    CNOlistTmp@inhibitors <- as.matrix(CNOlistTmp@inhibitors[, unique(colnames(CNOlistTmp@inhibitors))])
    NEMlistTmp$B <- NEMlistTmp$B[, sort(c(grepCtrls, grepStimuli, grep(paste(vertices, collapse = "|"), colnames(NEMlist$B))))]
    NEMlistTmp$pvals <- NEMlistTmp$pvals[, sort(c(grepCtrls, grepStimuli, grep(paste(vertices, collapse = "|"), colnames(NEMlist$B))))]
    #NEMlistTmp$signalStates <- NEMlist$signalStates[sort(c(grepCtrls, grepStimuli, grep(paste(vertices, collapse = "|"), rownames(NEMlist$signalStates)))), ]
    grepFcs <- colnames(computeFc(CNOlistTmp, NEMlistTmp$exprs))
    NEMlistTmp$fc <- NEMlist$fc[, which(colnames(NEMlist$fc) %in% grepFcs)]
    ## reduction of CNOlist:
    #CNOlistTmp <- cnoFromData(NEMlistTmp$exprs, stimuli, inhibitors, signals = pGraph$vertices[[i]])
    CNOlistTmp <- CNOlist
    CNOlistTmp@cues <- CNOlist@cues[, which(colnames(CNOlist@cues) %in% c(stimuli, inhibitors))]
    getRows <- which(rownames(CNOlistTmp@cues) %in% colnames(NEMlistTmp$exprs))
    CNOlistTmp@cues <- CNOlistTmp@cues[getRows, ]
    if (length(inhibitors) > 1) {
      CNOlistTmp@inhibitors <- CNOlist@inhibitors[getRows, which(colnames(CNOlist@inhibitors) %in% inhibitors)]
    } else {
      CNOlistTmp@inhibitors <- as.matrix(CNOlist@inhibitors[getRows, which(colnames(CNOlist@inhibitors) %in% inhibitors)])
      colnames(CNOlistTmp@inhibitors) <- inhibitors
    }
    CNOlistTmp@stimuli <- CNOlist@stimuli[, which(colnames(CNOlist@stimuli) %in% stimuli)]
    CNOlistTmp@stimuli <- CNOlistTmp@stimuli[getRows, ]
    CNOlistTmp@signals[[1]]<- CNOlist@signals[[1]][, which(colnames(CNOlist@signals[[1]]) %in% pGraph$vertices[[i]])]
    CNOlistTmp@signals[[1]] <- CNOlistTmp@signals[[1]][getRows, ]
    CNOlistTmp@signals[[2]]<- CNOlist@signals[[2]][, which(colnames(CNOlist@signals[[2]]) %in% pGraph$vertices[[i]])]
    CNOlistTmp@signals[[2]] <- CNOlistTmp@signals[[2]][getRows, ]
    ## optimize the subgraph
    #bString <- numeric(length(model$reacID))
    #bString[which(model$reacID %in% pGraph$edges[[i]])] <- 1
    #modelCut <- cutModel2(model, bString)

    PKNbot <- PKN
    PKNbot$reacID <- PKNbot$reacID[grep(paste("=", unlist(strsplit(pGraph$edges[[i]][1], "="))[2], sep = ""), PKNbot$reacID)]
    PKNbot$namesSpecies <- PKNbot$namesSpecies[which(PKNbot$namesSpecies %in% c(stimuli, inhibitors))]
    PKNbot$interMat <- PKNbot$interMat[grep(paste(c(stimuli,inhibitors), collapse = "|"), rownames(PKNbot$interMat)), grep(paste("=", unlist(strsplit(pGraph$edges[[i]], "="))[2], sep = ""), colnames(PKNbot$interMat))]
    PKNbot$notMat <- PKNbot$notMat[grep(paste(c(stimuli,inhibitors), collapse = "|"), rownames(PKNbot$notMat)), grep(paste("=", unlist(strsplit(pGraph$edges[[i]], "="))[2], sep = ""), colnames(PKNbot$notMat))]

    if (length(inhibitors) > 1) {
      makeSif(stimuli, inhibitors[-which(inhibitors %in% unlist(strsplit(pGraph$edges[[i]], "="))[2])], file = "temp.sif")
      PKNtop <- readSIF("temp.sif")
      unlink("temp.sif")
      PKNtmp <- PKNtop
      PKNtmp$reacID <- unique(c(PKNtop$reacID[grep(paste(stimuli, collapse = "|"), PKNtop$reacID)], PKNbot$reacID))
      PKNtmp$namesSpecies <- unique(c(PKNtop$namesSpecies, PKNbot$namesSpecies))
      PKNtmp$interMat <- cbind(rbind(PKNtop$interMat[, grep(paste(stimuli, collapse = "|"), colnames(PKNtop$interMat))], 0), PKNbot$interMat)
      PKNtmp$notMat <- cbind(rbind(PKNtop$notMat[, grep(paste(stimuli, collapse = "|"), colnames(PKNtop$notMat))], 0), PKNbot$notMat)
      
    } else {
      
      PKNtmp <- PKNbot
      
    }
    
    rownames(PKNtmp$interMat) <- rownames(PKNbot$interMat)
    rownames(PKNtmp$notMat) <- rownames(PKNbot$notMat)
    PKNtmp$interMat <- PKNtmp$interMat[, unique(colnames(PKNtmp$interMat))]
    PKNtmp$notMat <- PKNtmp$notMat[, unique(colnames(PKNtmp$notMat))]

    checkSignals(CNOlistTmp,PKNtmp)
    indices<-indexFinder(CNOlistTmp,PKNtmp,verbose=TRUE)
    NCNOindices<-findNONC(PKNtmp,indices,verbose=TRUE)
    NCNOcut<-cutNONC(PKNtmp,NCNOindices)
    indicesNCNOcut<-indexFinder(CNOlistTmp,NCNOcut)
    NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
    indicesNCNOcutComp<-indexFinder(CNOlistTmp,NCNOcutComp)
    NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate=100)
    resECNOlist<-residualError(CNOlistTmp)
    Fields4Sim<-prep4sim(NCNOcutCompExp)

    modelCut <- NCNOcutCompExp
    
    if (!is.null(initGates)) {
      initBstring <- numeric(length(modelCut$reacID))
      initBstring[which(modelCut$reacID %in% initGates)] <- 1
    }
    test <- 1
    if (!is.null(NEMlist$egenes) & test == 1) {
      relevant.sgenes <- colnames(CNOlistTmp@signals[[1]])
      relevant.egenes <- character()
      NEMlistTmp$egenes <- NULL
      for (sgene in relevant.sgenes) {
        relevant.egenes <- c(relevant.egenes, NEMlist$egenes[[sgene]])
        NEMlistTmp$egenes[[sgene]] <- NEMlist$egenes[[sgene]]
      }
      NEMlistTmp$exprs <- NEMlistTmp$exprs[relevant.egenes, ]
      NEMlistTmp$fc <- NEMlistTmp$fc[relevant.egenes, ]
      NEMlistTmp$B <- NEMlistTmp$B[relevant.egenes, ]
      NEMlistTmp$pvals <- NEMlistTmp$pvals[relevant.egenes, ]
    }
    
    result <- gaBinaryNemT1(CNOlistTmp,modelCut,initBstring=initBstring,sizeFac,NAFac,popSize,pMutation,maxTime,maxGens,stallGenMax,relTol,verbose,priorBitString,selPress,approach,NEMlistTmp,fit,targetBstring,elitism,inversion,graph,parameters,parallel,parallel2,selection)
    ## calculate the graph steady state for all experiments and the current graph (that also means experiments that were not used in the optimisation process)
    signalStates <- simulateStatesRecursive(CNOlist, modelCut, result$bString, NEMlist)
    #NEMlist$signalStates[, which(colnames(NEMlist$signalStates) %in% pGraph$vertices[[i]])] <- signalStates[, which(colnames(NEMlist$signalStates) %in% pGraph$vertices[[i]])]
    result$bString[-grep(paste("=", unlist(strsplit(pGraph$edges[[i]], "="))[2], sep = ""), modelCut$reacID)] <- 0
    graphs[[i]] <- modelCut$reacID[which(result$bString == 1)]
    
  }
  bString <- numeric(length(model$reacID))
  bString[which(model$reacID %in% unlist(graphs))] <- 1
  return(list(bString = bString, graphs = graphs, signalStates = NEMlist$signalStates))
}
