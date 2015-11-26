###### simple nem example:

source(".../cnopt.mod.R") # load all function

load(".../random.seed.RData") # if you want to recreate the example I had

library(CellNOptR)

stimuli <- "dummy"

while(length(stimuli) == 1) {

  dnf <- randomDnf(10, max.edges = 25, max.edge.size = 1, dag = TRUE)

  cues <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))))

  inputs <- unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)), "=")))

  outputs <- unique(gsub(".*=", "", dnf))

  stimuli <- inputs[which(!(inputs %in% outputs))]

}

inhibitors <- unique(c(inputs, outputs))
inhibitors <- inhibitors[-which(inhibitors %in% stimuli)]

plotDnf(dnf) # diamond heads denote actiavtion and negation in one arrow

sifMatrix <- NULL

for (i in dnf) {
  inputs2 <- unique(unlist(strsplit(gsub("=.*", "", i), "=")))
  output <- unique(gsub(".*=", "", i))
  for (j in inputs2) {
    j2 <- gsub("!", "", j)
    if (j %in% j2) {
      sifMatrix <- rbind(sifMatrix, c(j, 1, output))
    } else {
      sifMatrix <- rbind(sifMatrix, c(j2, -1, output))
    }
  }
}
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1, signals = NULL)

checkSignals(CNOlist,PKN)
indices<-indexFinder(CNOlist,PKN,verbose=TRUE)
NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
model<-expandNEM(NCNOcutComp, maxInputsPerGate=100)

bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

steadyState <- simulateStatesRecursive(CNOlist, model, bString)

while(any(apply(steadyState, 2, sd) == 0)) {

  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

  steadyState <- simulateStatesRecursive(CNOlist, model, bString)

}

plotDnf(model$reacID[as.logical(bString)])

NEMlist <- list()

NEMlist$exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 5)]

ERS <- computeFc(CNOlist, t(steadyState))

stimuli.pairs <- apply(apply(expand.grid(stimuli, stimuli), c(1,2), as.character), 1, paste, collapse = "_")

ERS <- ERS[, grep(paste(c(paste("Ctrl_vs", c(stimuli, inhibitors), sep = ""), paste(stimuli, "_vs_", stimuli, "_", rep(inhibitors, each = length(stimuli)), sep = ""), paste(stimuli.pairs, "_vs_", stimuli.pairs, "_", rep(inhibitors, each = length(stimuli.pairs)), sep = "")), collapse = "|"), colnames(ERS))]

NEMlist$fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)]
NEMlist$fc <- NEMlist$fc + rnorm(length(NEMlist$fc))

draw <- verbose <- T

initBstring <- reduceGraph(rep(0, length(model$reacID)), model) # start with empty graph; change 0 to 1 to start with PKN and remove reduce graph to start from fully connected (careful! might take long)

parameters <- list(cutOffs = c(0,1,0), scoring = c(0.13,0.05,0.25)) # can be used to calculate different score, but then method = NULL is required
sizeFac <- 10^-10 # zeta size penalty
method <- "s" # spearman rank correlation

parallel <- list(8, "rhskl9") # list(c(1,2,3,4), c("machine1", "machine2", "machine3", "machine4"))

res <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist, model=model, parallel=2, parallel2=1, seeds=1, initSeed = initBstring, sizeFac = sizeFac, method = method, verbose = verbose, draw = draw, max.steps = Inf, parameters = parameters)

par(mfrow=c(1,2), main = "ground truth network (left) and learned network (right)")

plotDnf(model$reacID[as.logical(bString)])

plotDnf(model$reacID[as.logical(res$bStrings[1, ])])

ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, res$bStrings[1, ])))

ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]

sum(ERS.res == ERS)/length(ERS) # accuracy of expected response scheme from learned network
