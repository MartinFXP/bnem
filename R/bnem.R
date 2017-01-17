#' Boolean Nested Effects Model main function
#' @param search Type of search heuristic. Either "greedy", "genetic" or "exhaustive". "greedy" uses a greedy algorithm to move through the local neighbourhood of a initial hyper-graph. "genetic" uses a genetic algorithm. "exhaustive" searches through the complete search space and is not recommended.
#' @param fc Foldchanges of gene expression values or equivalent input (normalized pvalues, logodds, ...). If left NULL, the gene expression data is used to calculate naive foldchanges.
#' @param exprs Optional normalized gene expression data.
#' @param pkn Prior knowledge network.
#' @param design Optional design matrix for the gene expression values, if available. If kept NULL, bnem needs either stimuli, inhibitors or a CNOlist object.
#' @param stimuli Character vector of stimuli names.
#' @param inhibitors Character vector of inhibitors.
#' @param signals Optional character vector of signals. Signals are S-genes, which can directly regulate E-genes. If left NULL, alls stimuli and inhibitors are defined as signals.
#' @param CNOlist CNOlist object if available.
#' @param model Model object including the search space, if available.
#' @param sizeFac Size factor penelizing the hyper-graph size.
#' @param NAFac factor penelizing NAs in the data.
#' @param NEMlist NEMlist object including the contrasts.
#' @param parameters atm not used
#' @param parallel Parallelize the search. An integer value specifies the number of threads on the local machine. A list object as in list(c(1,2,3), c("machine1", "machine2", "machine3")) specifies the threads distributed on different machines (local or others).
#' @param method Scoring method can be a correlation or distance measure. See ?cor and ?dist for details.
#' @param relFit atm not used
#' @param verbose TRUE gives additional information during the search.
#' @param reduce atm not used
#' @param initBstring Binary string of the initial hyper-graph.
#' @param popSize Popultaion size (only "genetic").
#' @param pMutation Probability for mutation (only "genetic").
#' @param maxTime Define a maximal time for the search.
#' @param maxGens Maximal number of generations (only "genetic").
#' @param stallGenMax Maximum number of stall generations (only "genetic").
#' @param relTol Score tolerance for networks defined as optimal but with a lower score as the real optimum (only "genetic").
#' @param priorBitString Binary string defining hyper-edges which are added to every hyper-graph. E.g. if you know hyper-edge 55 is definitly there and want to fix that set priorBitString[55] <- 1 (only "genetic").
#' @param selPress Selection pressure for the stochastic universal sampling (only "genetic").
#' @param approach atm not used
#' @param fit atm not used
#' @param targetBstring atmnot used
#' @param elitism Number of best hyper-graph transferred to the next generation  (only "genetic").
#' @param inversion Number of worst hyper-graphs for which their binary strings are inversed  (only "genetic").
#' @param parallel2 atm not used
#' @param selection "t" for tournament selection and "s" for stochastic universal sampling (only "genetic").
#' @param type atm not used
#' @param exhaustive If TRUE an exhaustive search is conducted if the genetic algorithm would take longer (only "genetic").
#' @param delcyc If TRUE deleted cycles in all hyper-graphs (recommended).
#' @param seeds atm not used
#' @param maxSteps Maximal number of steps (only "greedy").
#' @param node atm not used
#' @param absorpII Use a thrid absorption law.
#' @param draw If TRUE draws the network evolution.
#' @param prior Binary vector. A 1 specifies hyper-edges which should not be optimized (only "greedy").
#' @param ... additional low level parameters
#' @author Martin Pirkl
#' @seealso nem
#' @export
#' @import
#' CellNOptR
#' nem
#' snowfall
#' latticeExtra
#' @return List object including the optimized hyper-graph and its corresponding binary string.
#' @examples
#' library(bnem)
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signals = c("A", "B","C","D"))
#' model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnem(search = "greedy", model = model, CNOlist = CNOlist, fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"), parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
bnem <-
     function(search = "greedy",
 
              fc=NULL,
              exprs=NULL,
              egenes=NULL,
              pkn=NULL,
              design=NULL,
              stimuli=NULL,
              inhibitors=NULL,
              signals=NULL,
              CNOlist=NULL,
              model=NULL,
              sizeFac=10^-10,
              NAFac=1,
              parameters=list(cutOffs=c(0,1,0), scoring=c(0.1,0.2,0.9)),
              parallel = NULL,
              method = "s",
              approach = "fc",
              relFit = F,
              verbose = TRUE,
              reduce = TRUE,
             parallel2 = 1,

             initBstring = NULL,
             popSize = 100,
             pMutation = 0.5,
             maxTime = Inf,
             maxGens = Inf,
             stallGenMax = 10,
             relTol = 0.01,
             priorBitString = NULL,
             selPress = c(1.2,0.0001),
             fit = "linear",
             targetBstring = "none",
             elitism = NULL,
             inversion = NULL,
             selection = c("t"),
             type = "SOCK",
             exhaustive = FALSE,
             delcyc = TRUE,

             seeds = 1,
             maxSteps = Inf,
             node = NULL,
             absorpII = TRUE,
             draw = TRUE,
             prior = NULL,

             ...
             ) {
        if (is.null(model) | is.null(CNOlist)) {
            tmp <- preprocessInput(stimuli=stimuli,inhibitors=inhibitors,signals=signals,design=design,exprs=exprs,fc=fc,pkn=pkn,maxInputsPerGate=maxInputsPerGate)

            CNOlist <- tmp$CNOlist
            NEMlist <- tmp$NEMlist
            model <- tmp$model
        } else {
            NEMlist <- list()
            NEMlist$fc <- fc
            NEMlist$exprs <- exprs
            NEMlist$egenes <- egenes
        }

        if (search %in% c("greedy", "genetic", "exhaustive")) {

            if (search %in% "greedy") {
		res <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist, model=model, approach = approach, initSeed = initBstring, seeds = seeds, parameters = parameters, sizeFac = sizeFac, NAFac = NAFac, relTol = relTol, verbose = verbose, parallel=parallel, parallel2 = parallel2, relFit = relFit, method = method, max.steps = maxSteps, max.time = maxTime, node = node, absorpII = absorpII, draw = draw, prior = prior, ...)
		minSeed <- 1
		if (length(res$scores) > 1) {
                    for (i in 1:length(res$scores)) {
                        if (min(res$scores[[i]]) < min(res$scores[[minSeed]])) {
                            minSeed <- i
                        }
                    }
		}
		bString <- res$bStrings[minSeed, ]
		result <- list(graph = model$reacID[as.logical(bString)], bString = bString, bStrings = res$bStrings, scores = res$scores)
            }
            if (search %in% "genetic") {
		res <- gaBinaryNemT1(CNOlist=CNOlist, model=model,initBstring = initBstring,sizeFac = sizeFac, NAFac = NAFac,popSize = popSize,pMutation = pMutation,maxTime = maxTime,maxGens = maxGens, stallGenMax = stallGenMax,relTol = relTol,verbose = verbose,priorBitString = priorBitString,selPress = selPress,approach = approach,NEMlist=NEMlist,fit = fit,targetBstring = targetBstring,elitism = elitism,inversion = inversion,graph = draw,parameters = parameters,parallel = parallel,parallel2 = parallel2,selection = selection,relFit = relFit,method = method,type = type,exhaustive = exhaustive,delcyc = delcyc, ...)
		result <- list(graph = model$reacID[as.logical(res$bString)], bString = res$bString, bStrings = res$stringsTol, scores = res$stringsTolScores)
            }
            if (search %in% "exhaustive") {
		res <- exSearch(CNOlist=CNOlist,model=model,sizeFac=sizeFac,NAFac=NAFac,NEMlist=NEMlist,parameters=parameters, parallel = parallel, method = method, relFit = relFit, verbose = verbose, reduce = reduce, approach = approach, ...)
		result <- list(graph = model$reacID[as.logical(res$bString)], bString = res$bString, bStrings = res$bStrings, scores = res$scores)
            }
            return(result)
        } else {
            return("search must be be one of greedy, genetic or exhaustive")
        }
    }
