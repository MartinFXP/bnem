computeScoreNemT1 <- function(CNOlist,
                              model,
                              bString,
                              sizeFac = 10^-10,
                              NAFac = 1,
                              parameters = list(cutOffs = c(0,1,0), scoring = c(0.25,0.5,2)),
                              approach = "fc",
                              NEMlist = NULL,
                              tellme = 0,
                              sim = 0,
                              relFit = FALSE,
                              method = "s",
                              verbose = TRUE,
                              opt = "min"
                              ) {
  cutModel2 <- function (model, bString) {
    if (sum(bString == 1) > 0) {
      bs = as.logical(bString)
      newmodel <- list()
      if (is.null(dim(model$interMat))) {
        newmodel$interMat <- model$interMat[bs]
        newmodel$notMat <- model$notMat[bs]
        newmodel$reacID <- model$reacID[bs]
        newmodel$namesSpecies <- model$namesSpecies
      } else {
        newmodel$interMat <- model$interMat[, bs]
        newmodel$notMat <- model$notMat[, bs]
        newmodel$reacID <- model$reacID[bs]
        newmodel$namesSpecies <- model$namesSpecies
      }
    } else {
      newmodel <- model
    }
    return(newmodel)
  }

  CNOlist <- checkCNOlist(CNOlist)
  method <- checkMethod(method)
  ## if (is.null(simList) == TRUE) {
  ##   simList = prep4sim(model)
  ## }
  ## if (is.null(indexList) == TRUE) {
  ##   indexList = indexFinder(CNOlist, model)
  ## }
  if (sim == 0) {
    timeIndex = 2
    modelCut = cutModel2(model, bString)
    if (all(bString == 0)) {
      simResults <- simulateStatesRecursive(CNOlist = CNOlist, model = model, bString = bString, NEMlist)
    } else {
      simResults <- simulateStatesRecursive(CNOlist = CNOlist, model = modelCut, bString = (numeric(length(modelCut$reacID)) + 1), NEMlist)
    }
  }
  ## if (sim == 1) {
  ##   print("test2")
  ##   simResults <- simulateStatesRecursiveAdd(CNOlist = CNOlist, model = modelCut, bString = (numeric(length(modelCut$reacID)) + 1), NEMlist)
  ## }
  ## if (sim == 2) {
  ##   simListCut <- cutSimList(simList, bString)
  ##   simResults <- simulatorT1(CNOlist = CNOlist, model = modelCut, simList = simListCut, indexList = indexList)
  ## }
  ## nInTot <- length(which(model$interMat == -1))  does that work?
  nInTot <- length(unlist(strsplit(model$reacID, "\\+")))
  nTmp <- unlist(strsplit(model$reacID[which(bString == 1)], "\\+"))
  nInputsNeg <- length(grep("!", nTmp))
  nInputs <- length(nTmp) - nInputsNeg
  ## nInTot <- 1 # use this for no normalization
  sizePen <- sum((c(nInputs, nInputsNeg)/nInTot)*sizeFac)#/nrow(NEMlist$fc)
  Score <- getNemFit(simResults = simResults, CNOlist = CNOlist, model = modelCut, timePoint = timeIndex, NAFac = NAFac, sizePen = sizePen, simResultsT0 = NA, approach = approach, NEMlist = NEMlist, tellme = tellme, parameters = parameters, sim = sim, relFit = relFit, method = method, verbose = verbose, opt = opt)
  if ((class(CNOlist) == "CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  return(Score)
}
