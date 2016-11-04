gaRecon <- function (CNOlist, model, initBstring = NULL, # initBstring = TRUE
                           sizeFac = 1e-04, 
                           NAFac = 1, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500, 
                           stallGenMax = 100, elitism = NULL, relTol = 0.1, 
                           verbose = TRUE, priorBitString = NULL,
                           selPress = c(1.2,0.0001), # 1.2
                           approach = "fc",
                           NEMlist,
                           fit = "linear",
                           targetBstring = "none",
                           inversion = NULL,
                           graph = FALSE,
                           parameters = list(cutOffs = c(0.5,0.5), scoring = c(1,0.75), fitBest = 0),
                           parallel = NULL # parallel = N with N number of cores to use
                           ) {
  ## get all triples:
  triples <- list()
  stimuli <- colnames(CNOlist@stimuli)
  memory <- character()
  count <- 0
  for (i in colnames(CNOlist@inhibitors)) {
    for (j in colnames(CNOlist@inhibitors)) {
      if (i %in% j) { next() }
      for (k in colnames(CNOlist@inhibitors)) {
        if (i %in% k | j %in% k) { next() }
        if (paste(sort(c(i, j, k)), collapse = "_") %in% memory) { next() }
        count <- count + 1
        memory <- c(memory, paste(sort(c(i, j, k)), collapse = "_"))
        triples[[count]] <- sort(c(i, j, k))
      }
    }
  }
  triplesGates <- character()
  CNOlistOrg <- CNOlist
  NEMlistOrg <- NEMlist
  if (!is.null(parallel)) {
    require(snowfall)
    #maxCores <- parallel::detectCores()
    #sfSetMaxCPUs(number=parallel)
    if (is.list(parallel)) {
      if (length(parallel[[1]]) != length(parallel[[2]])) { stop("The nodes (second list object in parallel) and the number of cores used on every node (first list object in parallel) must be the same.") }
      hosts <- character()
      for (i in 1:length(parallel[[1]])) {
        hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
      }
      hosts <- as.list(hosts)
      sfInit(parallel=TRUE, socketHosts=hosts)
    } else {
      sfInit(parallel=TRUE, cpus=parallel)
    }
    sfExportAll()
    sfLibrary(CellNOptR)
    sfLibrary(hash)
  }
  cnolistStimuli <- apply(CNOlist@stimuli, 1, sum)
  cnolistInhibitors <- apply(CNOlist@inhibitors, 1, sum)
  cnolistCues <- apply(CNOlist@cues, 1, sum)
  maxStim <- max(cnolistStimuli)
  maxKd <- max(cnolistInhibitors)
  grepCtrl <- which(cnolistCues == 0)[1]
  grepStims <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors == 0))
  grepKds <- intersect(which(cnolistStimuli == 0), which(cnolistInhibitors != 0))
  grepStimsKds <- intersect(which(cnolistStimuli != 0), which(cnolistInhibitors != 0))
  for (i in 1:length(triples)) {
    inhibitors <- triples[[i]]
    CNOlist@cues <- CNOlistOrg@cues[, c(stimuli, inhibitors)]
    CNOlist@stimuli <- CNOlistOrg@stimuli[, stimuli]
    CNOlist@inhibitors <- CNOlistOrg@inhibitors[, inhibitors]
    CNOlist@signals[[1]] <- CNOlistOrg@signals[[1]][, c(inhibitors, "D1", "D2")]
    CNOlist@signals[[2]] <- CNOlistOrg@signals[[2]][, c(inhibitors, "D1", "D2")]
    cnolistStimuliTmp <- apply(CNOlist@stimuli, 1, sum)
    cnolistInhibitorsTmp <- apply(CNOlist@inhibitors, 1, sum)
    grepKdsTmp <- intersect(which(cnolistStimuliTmp == 0), which(cnolistInhibitorsTmp != 0))
    grepStimsKdsTmp <- intersect(which(cnolistStimuliTmp != 0), which(cnolistInhibitorsTmp != 0))
    # first make the "pkn":
    makeSif(stimuli, inhibitors, file = "PKNtemp.sif")
    PKN <- readSIF(sifFile="PKNtemp.sif")
    unlink("PKNtemp.sif")
    checkSignals(CNOlist,PKN)
    indices<-indexFinder(CNOlist,PKN,verbose=TRUE)
    NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
    NCNOcut<-cutNONC(PKN,NCNOindices)
    indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
    NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
    indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
    NCNOcutCompExp<-expandGates(NCNOcutComp)
    resECNOlist<-residualError(CNOlist)
    Fields4Sim<-prep4sim(NCNOcutCompExp)
    
    #egenes <- egeneSel(CNOlist, NEMlist, type, parameters)
    #NEMlist$exprs <- NEMlistOrg$exprs[unique(unlist(egenes)), ]
    NEMlist$exprs <- NEMlistOrg$exprs[, c(grepCtrl, grepStims, intersect(grepKdsTmp, which(cnolistInhibitors <= 3)))]
    NEMlist$exprs <- NEMlistOrg$exprs[, c(grepCtrl, grepStims, grepKdsTmp)]

    signals <- inhibitors

    CNOlist <- cnoFromData(NEMlist$exprs, stimuli = stimuli, inhibitors = inhibitors, signals = signals)
    CNOlist <- addDumpnodes(CNOlist)

    NEMlist$fc <- computeFc(CNOlist, NEMlist$exprs)

    initBstring<-rep(0,length(NCNOcutCompExp$reacID))
    CNOresult <- gaBinaryNemT1(
                               CNOlist=CNOlist,
                               NEMlist = NEMlist,
                               model=NCNOcutCompExp,
                               initBstring=initBstring,
                               popSize = popSize,
                               pMutation = pMutation,
                               maxTime = Inf,
                               maxGens = Inf,
                               stallGenMax = stallGenMax,
                               selPress = selPress,
                               sizeFac = sizeFac,
                               elitism = elitism,
                               relTol = relTol,
                               verbose = verbose,
                               fit = fit,
                               approach = approach,
                               parameters = parameters,
                               parallel = parallel
                               )

    triplesGates <- c(triplesGates, NCNOcutCompExp$reacID[which(CNOresult$bString == 1)])
    
  }

  triplesGates <- unique(triplesGates)

  return(triplesGates)

}
    
  
