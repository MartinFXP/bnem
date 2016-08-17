bnem <-
function(search = "greedy",

CNOlist,
NEMlist,
model,
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
graph = TRUE,
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
		res <- gaBinaryNemT1(CNOlist=CNOlist, model=model,initBstring = initBstring,sizeFac = sizeFac, NAFac = NAFac,popSize = popSize,pMutation = pMutation,maxTime = maxTime,maxGens = maxGens, stallGenMax = stallGenMax,relTol = relTol,verbose = verbose,priorBitString = priorBitString,selPress = selPress,approach = approach,NEMlist=NEMlist,fit = fit,targetBstring = targetBstring,elitism = elitism,inversion = inversion,graph = graph,parameters = parameters,parallel = parallel,parallel2 = parallel2,selection = selection,relFit = relFit,method = method,type = type,exhaustive = exhaustive,delcyc = delcyc, ...)
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
