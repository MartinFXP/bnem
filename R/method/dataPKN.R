dataPKN <- function (CNOlist, sizeFac = 1e-04, 
                      NAFac = 1, popSize = 100, pMutation = 0.5, maxTime = Inf, maxGens = Inf, 
                      stallGenMax = 100, elitism = NULL, relTol = 0.01, 
                      verbose = TRUE, maxSizeHashTable = 5000,
                      selPress = c("change", "1.2", "0.5", "3", "1"),
                      approach = "fc",
                      NEMlist,
                      fit = "nonlinear",
                      targetBstring = "none",
                      immigration = NULL,
                      singleEdge = FALSE,
                      graph = FALSE,
                      parameters = c(0.75, 1),
                      egenePars = c(log2(1.3), 0.05, 0)
                           ) {
  # get all triples:
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
  for (i in 922:length(triples)) {
    inhibitors <- triples[[i]]
    CNOlist@cues <- CNOlistOrg@cues[, c(stimuli, inhibitors)]
    CNOlist@stimuli <- CNOlistOrg@stimuli[, stimuli]
    CNOlist@inhibitors <- CNOlistOrg@inhibitors[, inhibitors]
    CNOlist@signals[[1]] <- CNOlistOrg@signals[[1]][, c(inhibitors, "Dump1", "Dump2")]
    CNOlist@signals[[2]] <- CNOlistOrg@signals[[2]][, c(inhibitors, "Dump1", "Dump2")]
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
                               immigration = immigration,
                               relTol = relTol,
                               verbose = verbose,
                               fit = fit,
                               graph = graph,
                               approach = approach,
                               targetBstring = targetBstring,
                               singleEdge = singleEdge,
                               parameters = parameters
                               )

    triplesGates <- c(triplesGates, NCNOcutCompExp$reacID[which(CNOresult$bString == 1)])
    
  }

  triplesGates <- unique(triplesGates)

}
    
  
