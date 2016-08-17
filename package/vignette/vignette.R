

library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/', fig.align='center', fig.show='hold')
options(formatR.arrow=TRUE,width=90)



X11.options(type="Xlib")
## install.packages("bnem_1.0.tar.gz")
library(bnem)



set.seed(2579)
## alternative seed and also a great song:
## set.seed(9247)

## to get the while loop started which makes sure we get a PKN with exactly two stimuli (not necessary):
stimuli <- "dummy"

while(length(stimuli) != 2) {

## random Boolean graph without cycles, maximal 25 edges, maximal edge size of 1 (normal graph) and maximal 10 S-genes:
  dnf <- randomDnf(10, max.edges = 25, max.edge.size = 1, dag = T) 
  
## all S-genes:
  cues <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))))

## parents:
  inputs <- unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)), "=")))

## children:
  outputs <- unique(gsub(".*=", "", dnf))

## parents which are no children are stimuli:
  stimuli <- inputs[which(!(inputs %in% outputs))]

}

inhibitors <- unique(c(inputs, outputs))
## S-genes which are no stimuli are inhibitors
inhibitors <- inhibitors[-which(inhibitors %in% stimuli)]



plotDnf(dnf, stimuli = stimuli)



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
## unlink("temp.sif")

## create dummy metainformation:
CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1, signals = NULL)

## extend the model:
model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100) 



## define a bitstring denoting the present and absent edges:
bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

## simulate S-gene states for all possible conditions:
steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)

## we find constitutively active S-genes with the folloing:
steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] <- steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] + CNOlist@inhibitors

## this while loop makes sure we get a gtn which actually affects all vertices and no vertices are constitutively active (not necessary):
while(any(apply(steadyState, 2, sd) == 0) | any(apply(steadyState2, 2, sd) == 0)) {

  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

  steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)

  steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] <- steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] + CNOlist@inhibitors

}



plotDnf(model$reacID[as.logical(bString)], stimuli = stimuli)



## this objects holds all data:
NEMlist <- list() 

## the expression data with 10 E-genes for each S-gene and 3 replicates (not used):
NEMlist$exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 3)]

## we calculate the foldchanges or expected S-gene response scheme (ERS) between certain condtion (e.g. control vs stimulation):
ERS <- computeFc(CNOlist, t(steadyState))

# the next step reduces the ERS to a sensible set of comparisons; e.g. we do not want to compare stimuli vs inhibition, but stimuli vs (stimuli,inhibition):
stimuli.pairs <- apply(apply(expand.grid(stimuli, stimuli), c(1,2), as.character), 1, paste, collapse = "_")

## this is the usual setup, but "computeFc" calculates a lot more contrasts, which can also be used if preferred:
ERS <- ERS[, grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""), paste(stimuli, "_vs_", stimuli, "_", rep(inhibitors, each = length(stimuli)), sep = ""), paste(stimuli.pairs, "_vs_", stimuli.pairs, "_", rep(inhibitors, each = length(stimuli.pairs)), sep = "")), collapse = "|"), colnames(ERS))]

## same as before with the expression values, we have 10 E-genes each and 3 replicates:
NEMlist$fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)]

## we add some Gaussian noise:
NEMlist$fc <- NEMlist$fc + rnorm(length(NEMlist$fc), 0, 1)

## some E-genes are negative regulated, hence we flip their foldchanges:
flip <- sample(1:nrow(NEMlist$fc), floor(0.33*row(NEMlist$fc)))
NEMlist$fc[flip, ] <- NEMlist$fc[flip, ]*(-1)

rownames(NEMlist$fc) <- paste(rownames(NEMlist$fc), 1:nrow(NEMlist$fc), sep = "_")
print(NEMlist$fc[1:6, c(1:3,(ncol(NEMlist$fc)-2):ncol(NEMlist$fc))])



## start with empty graph:
initBstring <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)
## initBstring <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)

## paralellize for several threads on one machine or multiple machines. See package "snowfall" for details.
parallel <- 2 # list(c(4,16,8,2), c("machine1", "machine2", "machine3", "machine4"))

## greedy search:
greedy <- bnem(search = "greedy",
            CNOlist=CNOlist,
            NEMlist=NEMlist,
            model=model,
            parallel=parallel,
            initBstring=initBstring,
            draw = FALSE,
            verbose = FALSE,
            maxSteps = Inf
            )

resString <- greedy$bString



par(mfrow=c(1,2))

## GTN:
plotDnf(model$reacID[as.logical(bString)], main = "GTN", stimuli = stimuli)

## optimum found:
plotDnf(model$reacID[as.logical(resString)], main = "greedy optimum", stimuli = stimuli)



## hyper-edge sensitivity and specificity:
print(sum(bString == 1 & resString == 1)/(sum(bString == 1 & resString == 1) + sum(bString == 1 & resString == 0)))
print(sum(bString == 0 & resString == 0)/(sum(bString == 0 & resString == 0) + sum(bString == 0 & resString == 1)))

## accuracy of expected response scheme from learned network should be high even though the network can look different:
ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, resString)))
ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]
print(sum(ERS.res == ERS)/length(ERS))



fitinfo <- validateGraph(CNOlist, NEMlist, model = model, bString = resString, Sgene = 1, Egenes = 1000, cexRow = 0.8, cexCol = 0.7, xrot = 45, disc = 0, Colv = T, Rowv = T, dendrogram = "both", bordercol = "lightgrey", aspect = "iso", sub = "")



## genetic algorithm:
genetic <- bnem(search = "genetic",
           parallel = parallel,
           CNOlist=CNOlist,
           NEMlist = NEMlist,
           model=model,
           initBstring=initBstring,
           popSize = 10,
           stallGenMax = 10,
           graph = FALSE,
           verbose = FALSE
           )

resString <- genetic$bString



## ## exhaustive search:
## exhaustive <- bnem(search = "exhaustive",
##            parallel = parallel,
##            CNOlist=CNOlist,
##            NEMlist = NEMlist,
##            model=model
##            )

## resString <- exRun$bString



## we do not force a DAG but do not allow repression:
dnf <- randomDnf(10, max.edges = 25, max.edge.size = 1, dag = F, negation = F) 
cues <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))))
inputs <- unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)), "=")))
outputs <- unique(gsub(".*=", "", dnf))
stimuli <- c(inputs[which(!(inputs %in% outputs))], cues[sample(1:length(cues), 2)])
inhibitors <- cues



both <- stimuli[which(stimuli %in% inhibitors)]
for (i in both) {
  dnf <- gsub(i, paste(i, "inhibit", sep = ""), dnf)
  dnf <- c(dnf, paste(i, "stim=", i, "inhibit", sep = ""))
  stimuli <- gsub(i, paste(i, "stim", sep = ""), stimuli)
  inhibitors  <- gsub(i, paste(i, "inhibit", sep = ""), inhibitors)
}



plotDnf(dnf, stimuli = stimuli)



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
## unlink("temp.sif")
CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1, signals = NULL)
model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100) 



bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)
steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)
steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] <- steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] + CNOlist@inhibitors
while(any(apply(steadyState, 2, sd) == 0) | any(apply(steadyState2, 2, sd) == 0)) {
  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)
  steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)
  steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] <- steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] + CNOlist@inhibitors
}
## we make sure the stimulations work:
bString[grep("stim", model$reacID)] <- 1
bString <- absorption(bString,model)



plotDnf(model$reacID[as.logical(bString)], stimuli = stimuli)



NEMlist <- list() 
NEMlist$exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 3)]
ERS <- computeFc(CNOlist, t(steadyState))
stimuli.pairs <- apply(apply(expand.grid(stimuli, stimuli), c(1,2), as.character), 1, paste, collapse = "_")
ERS <- ERS[, grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""), paste(stimuli, "_vs_", stimuli, "_", rep(inhibitors, each = length(stimuli)), sep = ""), paste(stimuli.pairs, "_vs_", stimuli.pairs, "_", rep(inhibitors, each = length(stimuli.pairs)), sep = "")), collapse = "|"), colnames(ERS))]
NEMlist$fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)]
NEMlist$fc <- NEMlist$fc + rnorm(length(NEMlist$fc), 0, 1)
flip <- sample(1:nrow(NEMlist$fc), floor(0.33*row(NEMlist$fc)))
NEMlist$fc[flip, ] <- NEMlist$fc[flip, ]*(-1)
rownames(NEMlist$fc) <- paste(rownames(NEMlist$fc), 1:nrow(NEMlist$fc), sep = "_")
print(NEMlist$fc[1:6, c(1:3,(ncol(NEMlist$fc)-2):ncol(NEMlist$fc))])



initBstring <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)
greedy2 <- bnem(search = "greedy",
            CNOlist=CNOlist,
            NEMlist=NEMlist,
            model=model,
            parallel=parallel,
            initBstring=initBstring,
            draw = FALSE,
            verbose = FALSE,
            maxSteps = Inf
            )
resString2 <- greedy2$bString



par(mfrow=c(1,2))
plotDnf(model$reacID[as.logical(bString)], main = "GTN", stimuli = stimuli)
plotDnf(model$reacID[as.logical(resString2)], main = "greedy optimum", stimuli = stimuli)



print(sum(bString == 1 & resString2 == 1)/(sum(bString == 1 & resString2 == 1) + sum(bString == 1 & resString2 == 0)))
print(sum(bString == 0 & resString2 == 0)/(sum(bString == 0 & resString2 == 0) + sum(bString == 0 & resString2 == 1)))
ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, resString2)))
ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]
print(sum(ERS.res == ERS)/length(ERS))



sessionInfo()


