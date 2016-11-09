## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/', fig.align='center', fig.show='hold')
options(formatR.arrow=TRUE,width=90)

## ----installandload---------------------------------------------------------------------
## install.packages("devtools")

library(devtools)

install_github("MartinFXP/B-NEM", quiet = T)

library(bnem)

## ----createpkn--------------------------------------------------------------------------
set.seed(2579)
## alternative seed and also a great song:
## set.seed(9247)

## to get the while loop started,
## which makes sure we get a PKN with exactly two stimuli (or what the amount you want):
stimuli <- "dummy"

while(length(stimuli) != 2) {

## random Boolean graph without cycles,
## maximal 25 edges, maximal edge size of 1 (normal graph) and maximal 10 S-genes:
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

## ----plotpkn, fig.width=4, fig.height=4, out.width='0.4\\linewidth'---------------------
plotDnf(dnf, stimuli = stimuli)

## ----pknextend--------------------------------------------------------------------------
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
write.table(sifMatrix, file = "temp.sif", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
## unlink("temp.sif")

## create metainformation (which S-genes are perturbed in which experiments):
CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors,
                        maxStim = 2, maxInhibit = 1, signals = NULL)

## extend the model:
model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100, verbose = T)

## ----definegtn--------------------------------------------------------------------------
## define a bitstring denoting the present and absent edges:
bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

## simulate S-gene states for all possible conditions:
steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)

## we find constitutively active S-genes with the following
## (they are boring, so we avoid them):
ind <- grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))
steadyState2[, ind] <- steadyState2[, ind] + CNOlist@inhibitors

## this while loop makes sure we get a gtn,
## which actually affects all vertices
## and no vertices are constitutively active:
while(any(apply(steadyState, 2, sd) == 0) | any(apply(steadyState2, 2, sd) == 0)) {

  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

  steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)

  ind <- grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))
  steadyState2[, ind] <- steadyState2[, ind] + CNOlist@inhibitors

}

## ----plotgtn, fig.width=4, fig.height=4, out.width='0.4\\linewidth'---------------------
plotDnf(model$reacID[as.logical(bString)], stimuli = stimuli)

## ----datasim----------------------------------------------------------------------------
## the expression data with 10 E-genes for each S-gene and 3 replicates:
exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 3)]

## we calculate the foldchanges (expected S-gene response scheme, ERS)
## between certain condtions (e.g. control vs stimulation):
ERS <- computeFc(CNOlist, t(steadyState))

# the next step reduces the ERS to a sensible set of comparisons;
## e.g. we do not want to compare stimuli vs inhibition,
## but stimuli vs (stimuli,inhibition):
stimcomb <- apply(expand.grid(stimuli, stimuli), c(1,2), as.character)
stimuli.pairs <- apply(stimcomb, 1, paste, collapse = "_")

## this is the usual setup,
## but "computeFc" calculates a lot more contrasts, which can also be used, if preferred:
ind <- grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""),
                    paste(stimuli, "_vs_", stimuli, "_",
                          rep(inhibitors, each = length(stimuli)), sep = ""),
                    paste(stimuli.pairs, "_vs_", stimuli.pairs, "_",
                          rep(inhibitors, each = length(stimuli.pairs)), sep = "")),
                  collapse = "|"), colnames(ERS))
ERS <- ERS[, ind]

## same as before with the expression values, we have 10 E-genes each and 3 replicates:
fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)]

## we add some Gaussian noise:
fc <- fc + rnorm(length(fc), 0, 1)

## in real applications some E-genes are negatively regulated,
## hence we flip some foldchanges:
flip <- sample(1:nrow(fc), floor(0.33*row(fc)))
fc[flip, ] <- fc[flip, ]*(-1)

## don't forget to set rownames (usually gene symbols, ensemble ids, entrez ids, ...)
rownames(fc) <- paste(rownames(fc), 1:nrow(fc), sep = "_")
print(fc[1:6, c(1:3,(ncol(fc)-2):ncol(fc))])

## ----localsearch------------------------------------------------------------------------
## we start with the empty graph:
initBstring <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)
## or a fully connected graph:
## initBstring <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)

## paralellize for several threads on one machine or multiple machines
## see package "snowfall" for details
parallel <- 2 # NULL for serialization
## or distribute to 30 threads on four different machines:
## parallel <- list(c(4,16,8,2), c("machine1", "machine2", "machine3", "machine4"))

## greedy search:
greedy <- bnem(search = "greedy",
            fc=fc,
            exprs=exprs, # not used, if fc is defined
            CNOlist=CNOlist,
            model=model,
            parallel=parallel,
            initBstring=initBstring,
            draw = FALSE, # TRUE: draw network at each step
            verbose = FALSE, # TRUE: print changed (hyper-)edges and score improvement
            maxSteps = Inf
            )

resString <- greedy$bString

## ----plotresult, fig.width=4, fig.height=4, out.width='0.4\\linewidth'------------------
par(mfrow=c(1,2))

## GTN:
plotDnf(model$reacID[as.logical(bString)], main = "GTN", stimuli = stimuli)

## greedy optimum:
plotDnf(model$reacID[as.logical(resString)], main = "greedy optimum", stimuli = stimuli)

## ----accuracy---------------------------------------------------------------------------
## hyper-edge sensitivity and specificity:
print(sum(bString == 1 & resString == 1)/
      (sum(bString == 1 & resString == 1) + sum(bString == 1 & resString == 0)))
print(sum(bString == 0 & resString == 0)/
      (sum(bString == 0 & resString == 0) + sum(bString == 0 & resString == 1)))

## accuracy of the expected response scheme (can be high even, if the networks differ):
ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, resString)))
ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]
print(sum(ERS.res == ERS)/length(ERS))

## ----visualfit, fig.width=20, fig.height=10, out.width='1.0\\linewidth'-----------------
fitinfo <- validateGraph(CNOlist, fc=fc, exprs=exprs, model = model, bString = resString,
                         Sgene = 1, Egenes = 1000, cexRow = 0.8, cexCol = 0.7, xrot = 45,
                         Colv = T, Rowv = T, dendrogram = "both", bordercol = "lightgrey",
                         aspect = "iso", sub = "")

## ----ga---------------------------------------------------------------------------------
## genetic algorithm:
genetic <- bnem(search = "genetic",
           fc=fc,
           exprs=exprs,
           parallel = parallel,
           CNOlist=CNOlist,
           model=model,
           initBstring=initBstring,
           popSize = 10,
           stallGenMax = 10,
           draw = FALSE,
           verbose = FALSE
           )

resString <- genetic$bString

## ----exhaustive-------------------------------------------------------------------------
## ## exhaustive search:
## exhaustive <- bnem(search = "exhaustive",
##            parallel = parallel,
##            CNOlist=CNOlist,
##            fc=fc,
##            exprs=exprs,
##            model=model
##            )

## resString <- exRun$bString

## ----createpkn2-------------------------------------------------------------------------
## we do not force a DAG but do not allow repression:
dnf <- randomDnf(10, max.edges = 25, max.edge.size = 1, dag = F, negation = F)
cues <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))))
inputs <- unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)), "=")))
outputs <- unique(gsub(".*=", "", dnf))
stimuli <- c(inputs[which(!(inputs %in% outputs))], cues[sample(1:length(cues), 2)])
inhibitors <- cues

## ----addstimuli-------------------------------------------------------------------------
both <- stimuli[which(stimuli %in% inhibitors)]
for (i in both) {
  dnf <- gsub(i, paste(i, "inhibit", sep = ""), dnf)
  dnf <- c(dnf, paste(i, "stim=", i, "inhibit", sep = ""))
  stimuli <- gsub(i, paste(i, "stim", sep = ""), stimuli)
  inhibitors  <- gsub(i, paste(i, "inhibit", sep = ""), inhibitors)
}

## ----plotpkn2, fig.width=6, fig.height=6, out.width='0.6\\linewidth'--------------------
plotDnf(dnf, stimuli = stimuli)

## ----pknextend2-------------------------------------------------------------------------
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
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
PKN <- readSIF("temp.sif")
## unlink("temp.sif")
CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors,
                        maxStim = 2, maxInhibit = 1, signals = NULL)
model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100, verbose = F)

## ----definegtn2-------------------------------------------------------------------------
bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)
steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)
ind <- grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))
steadyState2[, ind] <- steadyState2[, ind] + CNOlist@inhibitors
while(any(apply(steadyState, 2, sd) == 0) | any(apply(steadyState2, 2, sd) == 0)) {
  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)
  steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)
  ind <- grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))
  steadyState2[, ind] <- steadyState2[, ind] + CNOlist@inhibitors
}
## we make sure the stimulations work:
bString[grep("stim", model$reacID)] <- 1
bString <- absorption(bString,model)

## ----plotgtn2, fig.width=7, fig.height=6, out.width='0.6\\linewidth'--------------------
plotDnf(model$reacID[as.logical(bString)], stimuli = stimuli)

## ----datasim2---------------------------------------------------------------------------
exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 3)]
ERS <- computeFc(CNOlist, t(steadyState))
stmcomb <- apply(expand.grid(stimuli, stimuli), c(1,2), as.character)
stimuli.pairs <- apply(stimcomb, 1, paste, collapse = "_")
ind <- grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""),
                    paste(stimuli, "_vs_", stimuli, "_",
                          rep(inhibitors, each = length(stimuli)), sep = ""),
                    paste(stimuli.pairs, "_vs_", stimuli.pairs, "_",
                          rep(inhibitors, each = length(stimuli.pairs)), sep = "")),
                  collapse = "|"), colnames(ERS))
ERS <- ERS[, ind]
fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)]
fc <- fc + rnorm(length(fc), 0, 1)
flip <- sample(1:nrow(fc), floor(0.33*row(fc)))
fc[flip, ] <- fc[flip, ]*(-1)
rownames(fc) <- paste(rownames(fc), 1:nrow(fc), sep = "_")
print(fc[1:6, c(1:3,(ncol(fc)-2):ncol(fc))])

## ----localsearch2-----------------------------------------------------------------------
initBstring <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)
greedy2 <- bnem(search = "greedy",
            CNOlist=CNOlist,
            fc=fc,
            exprs=exprs,
            model=model,
            parallel=parallel,
            initBstring=initBstring,
            draw = FALSE,
            verbose = FALSE,
            maxSteps = Inf
            )
resString2 <- greedy2$bString

## ----plotresult2, fig.width=7, fig.height=6, out.width='0.6\\linewidth'-----------------
par(mfrow=c(1,2))
plotDnf(model$reacID[as.logical(bString)], main = "GTN", stimuli = stimuli)
plotDnf(model$reacID[as.logical(resString2)], main = "greedy optimum", stimuli = stimuli)

## ----accuracy2--------------------------------------------------------------------------
print(sum(bString == 1 & resString2 == 1)/
      (sum(bString == 1 & resString2 == 1) + sum(bString == 1 & resString2 == 0)))
print(sum(bString == 0 & resString2 == 0)/
      (sum(bString == 0 & resString2 == 0) + sum(bString == 0 & resString2 == 1)))
ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, resString2)))
ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]
print(sum(ERS.res == ERS)/length(ERS))

## ----definePos--------------------------------------------------------------------------
egenes <- list()

for (i in cues) {
    egenes[[i]] <- rownames(fc)[grep(i, rownames(fc))]
}

initBstring <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)
greedy2 <- bnem(search = "greedy",
            CNOlist=CNOlist,
            fc=fc,
            exprs=exprs,
            egenes=egenes,
            model=model,
            parallel=parallel,
            initBstring=initBstring,
            draw = FALSE,
            verbose = FALSE,
            maxSteps = Inf
            )
resString3 <- greedy2$bString

## ----accuracy3--------------------------------------------------------------------------
print(sum(bString == 1 & resString3 == 1)/
      (sum(bString == 1 & resString3 == 1) + sum(bString == 1 & resString3 == 0)))
print(sum(bString == 0 & resString3 == 0)/
      (sum(bString == 0 & resString3 == 0) + sum(bString == 0 & resString3 == 1)))
ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, resString3)))
ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]
print(sum(ERS.res == ERS)/length(ERS))

## ----residuals, fig.width=7, fig.height=5, out.width='0.9\\linewidth'-------------------
residuals <- findResiduals(resString3, CNOlist, model, fc, verbose = F) # verbose = TRUE plots the residuals matrices

## ----loadbcrdata------------------------------------------------------------------------
data(bcr)
head(fc)

## ----bcrpkn-----------------------------------------------------------------------------
negation <- F # what happens if we allow negation?
sifMatrix <- numeric()
for (i in "BCR") {
  sifMatrix <- rbind(sifMatrix, c(i, 1, c("Pi3k")))
  sifMatrix <- rbind(sifMatrix, c(i, 1, c("Tak1")))
  if (negation) {
    sifMatrix <- rbind(sifMatrix, c(i, -1, c("Pi3k")))
    sifMatrix <- rbind(sifMatrix, c(i, -1, c("Tak1")))
  }
}
for (i in c("Pi3k", "Tak1")) {
  for (j in c("Ikk2", "p38", "Jnk", "Erk", "Tak1", "Pi3k")) {
    if (i %in% j) { next() }
    sifMatrix <- rbind(sifMatrix, c(i, 1, j))
    if (negation) {
      sifMatrix <- rbind(sifMatrix, c(i, -1, j))
    }
  }
}
for (i in c("Ikk2", "p38", "Jnk")) {
  for (j in c("Ikk2", "p38", "Jnk")) {
    if (i %in% j) { next() }
    sifMatrix <- rbind(sifMatrix, c(i, 1, j))
    if (negation) {
      sifMatrix <- rbind(sifMatrix, c(i, -1, j))
    }
  }
}

write.table(sifMatrix, file = "temp.sif", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

## ----bcrmeta----------------------------------------------------------------------------
CNOlist <- dummyCNOlist(stimuli = "BCR",
                        inhibitors = c("Tak1", "Pi3k", "Ikk2", "Jnk", "p38", "Erk"),
                        maxStim = 1, maxInhibit = 3)

model <- preprocessing(CNOlist, PKN)

## ----bcra-------------------------------------------------------------------------------
initBstring <- rep(0, length(model$reacID))
ga <- bnem(search = "genetic",
               fc=fc,
               CNOlist=CNOlist,
               model=model,
               parallel=2,
               initBstring=initBstring,
               draw = FALSE,
               verbose = FALSE
               )
print(min(ga$scores))

## ----bcrgreedy--------------------------------------------------------------------------
initBstring <- rep(0, length(model$reacID))
greedy <- bnem(search = "greedy",
               fc=fc,
               CNOlist=CNOlist,
               model=model,
               parallel=2,
               initBstring=initBstring,
               draw = FALSE,
               verbose = FALSE
               )
print(min(greedy$scores[[1]]))

## ----plotresultbcr, fig.width=7, fig.height=5, out.width='0.6\\linewidth'---------------
par(mfrow=c(1,2))
plotDnf(PKN$reacID, main = "PKN", stimuli = "BCR")
plotDnf(ga$graph, main = "genetic optimum", stimuli = "BCR")
plotDnf(greedy$graph, main = "greedy optimum", stimuli = "BCR")

## ----sessioninfo------------------------------------------------------------------------
sessionInfo()

