###### this is a documented guide to allow anyone to do a B-NEM analysis on their own dataset:

X11.options(type="Xlib")

source(".../cnopt.mod.R") # load all functions

library(CellNOptR) # CNO package required

# first we define vertices which can be stimulated and vertices whcih can be inhibited. currently there is no overlap allowed. however this can be circumvented if for one vertice which can be both an additional stimuli vertice is introduced which can only go into the vertice which can be inhibited.

stimuli <- "dummy" # just to get the while loop started

while(length(stimuli) != 2) { # we want exactly two stimuli

  dnf <- randomDnf(10, max.edges = 100, max.edge.size = 1, dag = TRUE) # create random disjunctive normal form (dnf) of 10 literals (=vertices, S-genes) with at max 1 parent per edge and 25 edges and no cycles

  cues <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))) # get all the vertices/S-genes

  inputs <- unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)), "="))) # which vertices are parents

  outputs <- unique(gsub(".*=", "", dnf)) # which vertices are children

  stimuli <- inputs[which(!(inputs %in% outputs))] # vertices without parents are stimuli

}

inhibitors <- unique(c(inputs, outputs))
inhibitors <- inhibitors[-which(inhibitors %in% stimuli)] # we want all vertices not stimuli as inhibitors

par(mfrow=c(1,2))

plotDnf(dnf, legend = 1, stimuli = stimuli, inhibitors = inhibitors, signals = c(stimuli, inhibitors)) # diamond heads denote activation and negation in one arrow; use plotDnf(legend=1)

plotDnf(dnf, legend = 1, stimuli = stimuli, inhibitors = inhibitors, signals = c(stimuli, inhibitors), simulate = list(stimulated = stimuli[1], inhibited = inhibitors[1]))

# next we need to build the model from a prior network

sifMatrix <- NULL # we want to load a prior knowledge network (pkn) from a sif file as it is done in the CNO package

for (i in dnf) { # we take the dnf and write the edges in sif format:
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
PKN <- readSIF("temp.sif") # we load the pkn (load your own pkn sif file for your data)
unlink("temp.sif")

CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1, signals = NULL) # this is used for meta information, e.g. all possible conditions with at max two stimuli active and one vertice inhibited

checkSignals(CNOlist,PKN) # this is used as in the cno package to create the model with all possible hyper-edges:
indices<-indexFinder(CNOlist,PKN,verbose=TRUE)
NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
model<-expandNEM(NCNOcutComp, maxInputsPerGate=100)

plotDnf(model$reacID[-grep("\\+", model$reacID)])

bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model) # we define a random ground truth network (gtn)

steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString) # we simulate the steady states for all possible conditions

steadyState[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState))] <- steadyState[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState))] + CNOlist@inhibitors

while(any(apply(steadyState, 2, sd) == 0) | any(apply(steadyState2, 2, sd) == 0)) { # this while loop makes sure we get a gtn which actually affects all vertices and no vertices are constitutively active

  bString <- absorption(sample(c(0,1), length(model$reacID), replace = T), model)

  steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model, bString)

  steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] <- steadyState2[, grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))] + CNOlist@inhibitors

}

plotDnf(model$reacID[as.logical(bString)], stimuli = stimuli) # let's look at our gtn

NEMlist <- list() # this list includes the real data (expression values and foldchanges)

NEMlist$exprs <- t(steadyState)[rep(1:ncol(steadyState), 10), rep(1:nrow(steadyState), 5)] # every vertices directly affects the transcription of 10 mRNA products or E-genes

ERS <- computeFc(CNOlist, t(steadyState)) # we calculate the foldchanges or expected S-gene response scheme (ERS) between certain condtion (e.g. control vs stimulation)

stimuli.pairs <- apply(apply(expand.grid(stimuli, stimuli), c(1,2), as.character), 1, paste, collapse = "_") # the next step reduce the ERS to a sensible set of comparisons; e.g. we do not want to compare stimuli vs inhibition, but stimuli vs (stimuli,inhibition)

ERS <- ERS[, grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""), paste(stimuli, "_vs_", stimuli, "_", rep(inhibitors, each = length(stimuli)), sep = ""), paste(stimuli.pairs, "_vs_", stimuli.pairs, "_", rep(inhibitors, each = length(stimuli.pairs)), sep = "")), collapse = "|"), colnames(ERS))] # this is the usual setup. But "computeFc" calculates a lot mor contrasts, which can also be used if preferred.

NEMlist$fc <- ERS[rep(1:nrow(ERS), 10), rep(1:ncol(ERS), 3)] # same as before with the expression values, we have 10 E-genes each and 3 replicates
NEMlist$fc <- NEMlist$fc + rnorm(length(NEMlist$fc), 0, 1) # we add Gaussian noise
flip <- sample(1:nrow(NEMlist$fc), floor(0.33*nrow(NEMlist$fc)))
NEMlist$fc[flip, ] <- NEMlist$fc[flip, ]*(-1) # some E-genes are negatively regulated
rownames(NEMlist$fc) <- paste(rownames(NEMlist$fc), 1:nrow(NEMlist$fc), sep = "_")

# If you want to do a real world analysis, you can mostly do everything except the data generation as described above. Obviously you have to provide your own pkn and the data is not coming from a known network. The NEMlist$exprs/fc objects correspond to your data. The exprs data is not used (officially). The fc data are your foldchanges from limma/DEseq/edgeR/"favourite differential expression tool". It is important that you name your colnames correctly. E.g. let's assume one of your stimulations is called A. If you have a contrast "A - control", you have to give the column the name "Ctrl_vs_A". If B is another stimulation and C an inhibitor vertice the contrast "(B,C) - B" (the condition with B stimulated and C inhibited against only B stimulated) must be named "B_vs_B_C". For two stimulations you must have the name "A_B_vs_A_B_C". The alphabetical order is also very important. So "B_A_vs_B_A_C" is incorrect. Or if D is another inhibitor "A_vs_A_D_C" or "A_B_vs_A_B_D_C" are also incorrect. Look at colnames(NEMlist$fc) for examples, but have in mind what is a stimuli and what an inhibitor.

initBstring <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist) # start with empty graph; change 0 to 1 to start with PKN (usually takes long, but can get better results.) Use different (random) starting networks to reduce the chance of a local optimum. E.g. you can draw a random start and also use its opposite "1 - initBstring" to cover more space efficiently.

parallel <- 8 # list(c(4,16,8,2), c("machine1", "machine2", "machine3", "machine4")) # parallel calculation of edge improvements

# try greedy search:

res <- localSearch(
         CNOlist=CNOlist,
         NEMlist=NEMlist,
         model=model,
         parallel=parallel,
         initSeed=bString#initBstring
         )

par(mfrow=c(1,2), main = "ground truth network (left) and learned network (right)") # plot the result vs the gtn

plotDnf(model$reacID[as.logical(bString)])

plotDnf(model$reacID[as.logical(res$bStrings[1, ])])

sum(bString == 1 & res$bStrings[1, ] == 1)/(sum(bString == 1 & res$bStrings[1, ] == 1) + sum(bString == 1 & res$bStrings[1, ] == 0)) # hyper-edge sensitivity

sum(bString == 0 & res$bStrings[1, ] == 0)/(sum(bString == 0 & res$bStrings[1, ] == 0) + sum(bString == 0 & res$bStrings[1, ] == 1)) # hyper-edge specificity

# the algorithm looks for the smallest best fitting network; thus it can be that the GTN is actually not the best network having that ERS and there are smaller but equivalent networks

ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, res$bStrings[1, ])))

ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]

sum(ERS.res == ERS)/length(ERS) # accuracy of expected response scheme from learned network should be high even though the network can look different

# lets look at the data and how well it fits the resolved network:

geneLists <- list()
par(ask=T)
for (i in 1:ncol(CNOlist@signals[[1]])) {
  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = res$bStrings[1, ], Sgene = i, Egenes = 1000, parameters=parameters, cexRow = 0.8, soft = T, cexCol = 0.7, xrot = 45, method = method, disc = 0, Colv = T, Rowv = T, dendrogram = "both", bordercol = "grey", aspect = "iso", sub = "")
  dev.print("temp.pdf", device = pdf, width = 40, height = 10) # take a closer look
}; names(geneLists) <- colnames(CNOlist@signals[[1]]); par(ask=F)

# or genetic algorithm:

gaRun <- gaBinaryNemT1(
           parallel = parallel,
           CNOlist=CNOlist,
           NEMlist = NEMlist,
           model=model,
           initBstring=initBstring,
           popSize = 100,
           stallGenMax = 10
           )

par(mfrow=c(1,2), main = "ground truth network (left) and learned network (right)") # plot the result vs the gtn

plotDnf(model$reacID[as.logical(bString)])

plotDnf(model$reacID[as.logical(gaRun$bString)])

sum(bString == 1 & gaRun$bString == 1)/(sum(bString == 1 & gaRun$bString == 1) + sum(bString == 1 & gaRun$bString == 0)) # hyper-edge sensitivity

sum(bString == 0 & gaRun$bString == 0)/(sum(bString == 0 & gaRun$bString == 0) + sum(bString == 0 & gaRun$bString == 1)) # hyper-edge specificity

# the algorithm looks for the smallest best fitting network; thus it can be that the GTN is actually not the best network having that ERS and there are smaller but equivalent networks

ERS.res <- computeFc(CNOlist, t(simulateStatesRecursive(CNOlist, model, gaRun$bString)))

ERS.res <- ERS.res[, which(colnames(ERS.res) %in% colnames(ERS))]

sum(ERS.res == ERS)/length(ERS) # accuracy of expected response scheme from learned network should be high even though the network can look different

# lets look at the data and how well it fits the resolved network:

geneLists <- list()
par(ask=T)
for (i in 1:ncol(CNOlist@signals[[1]])) {
  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = gaRun$bString, Sgene = i, Egenes = 10, parameters=parameters, cexRow = 0.8, soft = T, cexCol = 0.7, xrot = 45, method = method, disc = 0, Colv = T, Rowv = T, dendrogram = "both", bordercol = "grey", aspect = "iso", sub = "")
  dev.print("temp.pdf", device = pdf, width = 20, height = 15) # take a closer look
}

names(geneLists) <- colnames(CNOlist@signals[[1]])

par(ask=F)

# one important parameters, which should be trained before the final optimization starts is zeta. This controls the sparseness of the solution and is set to 10^-10 by default. You can change that with "szeFac=a" with any number a in the optimization functions localSearch or gaBinaryNemT1.

# have fun with your own analysis

