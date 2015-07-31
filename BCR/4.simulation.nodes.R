##########################################################
###                                                    ###
###         number of nodes                            ###
###                                                    ###
##########################################################

######################## from here on in all seriousness
X11.options(type="Xlib")

library(CellNOptR)

source("Boutros10.svn/method/scripts/cnopt.mod.R")

############## first define the set of experiments (=data) available:

stimuli <- paste("S", 1:6, sep = "")
inhibitors <- c(paste("I", 11:16, sep = ""), paste("I", 21:26, sep = ""), paste("I", 31:36, sep = ""), paste("I", 41:46, sep = ""))
signals <- c(stimuli, inhibitors)

sifMatrix <- numeric()

for (i in stimuli) {
  for (j in inhibitors[1:6]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(i, "-1", j))
    } else {
      sifMatrix <- rbind(sifMatrix, c(i, "1", j))
    }
  }
}
for (i in inhibitors[1:6]) {
  for (j in inhibitors[7:12]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(i, "-1", j))
    } else {
      sifMatrix <- rbind(sifMatrix, c(i, "1", j))
    }
  }
}
for (j in inhibitors[7:12]) {
  for (k in inhibitors[13:18]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(j, "-1", k))
    } else {
      sifMatrix <- rbind(sifMatrix, c(j, "1", k))
    }
  }
}
for (k in inhibitors[13:18]) {
  for (l in inhibitors[19:24]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(k, "-1", l))
    } else {
      sifMatrix <- rbind(sifMatrix, c(k, "1", l))
    }
  }
}

write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

length(PKN$reacID)

## PKN <- readSIF("Kube011_BCR_CD40_Inhibitoren/publication/PKNjulio.sif")
## stimuli <- c("il6", "tgfa", "ins", "tnf", "il1")
## inhibitors <- unique(gsub("!", "", unlist(strsplit(PKN$reacID, "="))))
## inhibitors <- inhibitors[-which(inhibitors %in% stimuli)]

CNOlist <- dummyCNOlist(stimuli, inhibitors, maxStim = 2, maxInhibit = 1) # length(inhibitors))

checkSignals(CNOlist,PKN)
indices<-indexFinder(CNOlist,PKN,verbose=TRUE)
NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate=2)
resECNOlist<-residualError(CNOlist)
Fields4Sim<-prep4sim(NCNOcutCompExp)

model <- NCNOcutCompExp

bString <- numeric(length(model$reacID))
bString[-grep("\\+", model$reacID)] <- 1

edgecol <- rep("black", length(model$reacID[-grep("\\+", model$reacID)]))
edgecol[grep("!", model$reacID[-grep("\\+", model$reacID)])] <- "red"

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_super_pkn.pdf", height = 20, width = 10)

plotDnf(model$reacID[-grep("\\+", model$reacID)], legend = 0, width = 1, edgecol = edgecol, stimuli = stimuli)

dev.off()

length(model$reacID)

######################## do fundamental super pkn sim:

source("Boutros10.svn/method/scripts/cnopt.mod.R")

method <- "s"
parameters <- list(cutOffs = c(0,1,0), scoring = c(0.1,0.9,2))
sizeFac <- 10^-10

verbose <- TRUE
popSize <- 100
stallGenMax <- 10
maxGens <- Inf
noise <- 0.1

parallel <- 12#list(4, "rhskl4")

count <- 0

CNOinput <- list()

CNOresults <- list()

TN.list <- list()

for (n in rev(round(seq(10, 15, (30-10)/4)))) {
  ## n <- 30
  count <- count + 1
  CNOinput[[count]] <- list()
  CNOresults[[count]] <- list()
  for (c in 1:10) {
    ## c <- 1
    if (n == 30) {
      core <- c(paste("S", sample(1:6, 1), sep = ""), paste("I", 1:4, sample(1:6, 4, replace = T), sep = ""))
      core <- c(core, sample(signals[-which(signals %in% core)], n - 5))
      nodes <- core
      notnodes <- signals[-which(signals %in% core)]
      bString <- numeric(length(model$reacID))
      bString[grep(paste(nodes, collapse = "|"), model$reacID)] <- 1
      if (n < 30) {
        bString[grep(paste(notnodes, collapse = "|"), model$reacID)] <- 0
      }
      model1 <- cutModel(model, bString)
      TN <- numeric(length(model1$reacID))
      ## TN <- sample(c(0,1), length(model1$reacID), replace = T) # 50/50 chance
      ## TN[sample(1:length(model$reacID), 50)] <- 1 # or only size 50? but size 50 is larger than the pkn of a size 10 node pkn!!!
      TN[sample(1:length(model$reacID), ceiling(length(model$reacID)*0.1))] <- 1 # 10% network size
      if (length(grep(paste(paste("S", 1:6, sep = ""), collapse = "|"), model1$reacID[as.logical(TN)])) == 0) {
        TN[sample(grep(paste(stimuli[which(stimuli %in% nodes)], collapse = "|"), model1$reacID), sample(1:length(grep(paste(stimuli[which(stimuli %in% nodes)], collapse = "|"), model1$reacID)), 1))] <- 1
      }
      TN <- reduceGraph(TN, model1, CNOlist)
      TN.list[[c]] <- model1$reacID[as.logical(TN)]
    } else {
      TN.graph <- TN.list[[c]]
      nodes.tmp <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(TN.graph, "\\+")), "="))))
      core <- gsub("!", "", unlist(strsplit(unlist(strsplit(TN.graph[sample(grep("S", TN.graph), 1)], "=")), "\\+")))
      core <- c(core, nodes.tmp[sample(grep("I2", nodes.tmp), 1)], nodes.tmp[sample(grep("I3", nodes.tmp), 1)], nodes.tmp[sample(grep("I4", nodes.tmp), 1)])
      core <- c(core, sample(nodes.tmp[-which(nodes.tmp %in% core)], n - length(core)))
      nodes <- core
      notnodes <- signals[-which(signals %in% core)]
      bString <- numeric(length(model$reacID))
      bString[grep(paste(nodes, collapse = "|"), model$reacID)] <- 1
      if (n < 30) {
        bString[grep(paste(notnodes, collapse = "|"), model$reacID)] <- 0
      }
      model1 <- cutModel(model, bString)
      TN <- numeric(length(model1$reacID))
      TN[which(model1$reacID %in% TN.graph)] <- 1
    }
    TN.graph <- model1$reacID[as.logical(TN)]
    SimResults <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = TN)
    nodes <- sort(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(TN.graph, "=")), "\\+")))))
    notnodes <- signals[-which(signals %in% nodes)]
    SimResults <- SimResults[, grep(paste(nodes, collapse = "|"), colnames(SimResults))]
    SimResults <- SimResults[grep(paste(c("Ctrl", nodes), collapse = "|"), rownames(SimResults)), ]
    if (n < 30) {
      SimResults <- SimResults[-grep(paste(notnodes, collapse = "|"), rownames(SimResults)), ]
    }
    NEMlist <- list()
    NEMlist$exprs <- t(SimResults[rep(1:nrow(SimResults), 3), rep(1:ncol(SimResults), 10)])
    NEMlist$exprs <- NEMlist$exprs[, order(colnames(NEMlist$exprs))]
    colnames(NEMlist$exprs) <- paste(paste(sort(rep(rownames(SimResults), 3)), "_rep", 1:3, sep = ""), "_run1", sep = "")
    noisy <- sample(1:length(NEMlist$exprs), floor(noise*length(NEMlist$exprs)))
    NEMlist$exprs[noisy] <- 1 - NEMlist$exprs[noisy]
    NEMlist$fc <- computeFcII(NEMlist$exprs, stimuli[which(stimuli %in% nodes)], inhibitors[which(inhibitors %in% nodes)], paste("rep", 1:3, sep = ""), "run1")
    colnames(NEMlist$fc) <- gsub("_r.*", "", colnames(NEMlist$fc))
    relevant <- numeric()
    for (i in stimuli[which(stimuli %in% nodes)]) {
      relevant <- c(relevant, grep(paste("Ctrl_vs_", i, "$", sep = ""), colnames(NEMlist$fc)))
      relevant <- c(relevant, grep(paste("^", i, "_vs_", i, "_I[1-4][1-6]", sep = ""), colnames(NEMlist$fc)))
      for (j in stimuli[which(stimuli %in% nodes)]) {
        relevant <- c(relevant, grep(paste("Ctrl_vs_", paste(sort(c(i, j)), collapse = "_"), "$", sep = ""), colnames(NEMlist$fc)))
        relevant <- c(relevant, grep(paste("^", paste(sort(c(i, j)), collapse = "_"), "_vs_", paste(sort(c(i, j)), collapse = "_"), "_I[1-4][1-6]", sep = ""), colnames(NEMlist$fc)))
      }
    }
    for (i in inhibitors[which(inhibitors %in% nodes)]) {
      relevant <- c(relevant, grep(paste("Ctrl_vs_", i, "$", sep = ""), colnames(NEMlist$fc)))
    }
    relevant <- unique(relevant)
    if (length(relevant) > 0) {
      NEMlist$fc <- NEMlist$fc[, relevant]
    }
    initBstring <- reduceGraph(rep(0,length(model1$reacID)), model)
    CNOinput[[count]][[c]] <- list(model1 = model1, TN = TN)
    print(paste("size: ", n, " sample: ", c, sep = ""))
    dim(NEMlist$fc)
    start <- Sys.time()
    CNOresults[[count]][[c]] <- gaBinaryNemT1(parallel=parallel,CNOlist=CNOlist,NEMlist=NEMlist,model=model1,initBstring=initBstring,popSize = popSize, maxTime = Inf, maxGens = maxGens, stallGenMax = stallGenMax, elitism = ceiling(popSize*0.01),inversion = ceiling(popSize*0.01),verbose=TRUE,parameters = parameters,sizeFac = sizeFac,method = method)
    end <- Sys.time()
    print(end - start)
  }
}

source("temp.R")

computeScoreNemT1(CNOlist, model = model, TN, simList = NULL, indexList = NULL, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

par(mfrow=c(1,2))
plotDnf(model1$reacID[as.logical(TN)], CNOlist = CNOlist, legend = 0)
plotDnf(model1$reacID[as.logical(RN)], CNOlist = CNOlist, legend = 0)

## save(CNOlist, CNOresults, CNOinput, file = "superpkn_nodes3.RData")

sens <- spec <- time <- hyeds <- matrix(0, length(CNOresults), length(CNOresults[[1]]))

for (i in 1:length(CNOresults)) {
  for (j in 1:length(CNOresults[[i]])) {
    model1 <- CNOinput[[i]][[j]]$model1
    nodes <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(model1$reacID, "\\+")), "="))))
    hyeds[i, j] <- length(model1$reacID)
    TN <- CNOinput[[i]][[j]]$TN
    RN <- CNOresults[[i]][[j]]$bString
    ERS.RN <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = RN)
    ERS.TN <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = TN)
    FC.RN <- computeFc(CNOlist, t(ERS.RN))
    FC.TN <- computeFc(CNOlist, t(ERS.TN))
    relevant <- numeric()
    for (k in stimuli[which(stimuli %in% nodes)]) {
      relevant <- c(relevant, grep(paste("Ctrl_vs_", k, "$", sep = ""), colnames(FC.RN)))
      relevant <- c(relevant, grep(paste("^", k, "_vs_", k, "_I[0-9]", "$", sep = ""), colnames(FC.RN)))
      relevant <- c(relevant, grep(paste("^", k, "_vs_", k, "_I[1-4][0-9]", "$", sep = ""), colnames(FC.RN)))
      for (l in stimuli[which(stimuli %in% nodes)]) {
        relevant <- c(relevant, grep(paste("Ctrl_vs_", paste(sort(c(k, l)), collapse = "_"), "$", sep = ""), colnames(FC.RN)))
        relevant <- c(relevant, grep(paste("^", paste(sort(c(k, l)), collapse = "_"), "_vs_", paste(sort(c(k, l)), collapse = "_"), "_I[6-9]", sep = ""), colnames(FC.RN)))
        relevant <- c(relevant, grep(paste("^", paste(sort(c(k, l)), collapse = "_"), "_vs_", paste(sort(c(k, l)), collapse = "_"), "_I[1-4][0-9]", sep = ""), colnames(FC.RN)))
      }
    }
    for (k in inhibitors[which(inhibitors %in% nodes)]) {
      relevant <- c(relevant, grep(paste("Ctrl_vs_", k, "$", sep = ""), colnames(FC.RN)))
    }
    relevant <- unique(relevant)
    FC.TN <- FC.TN[, relevant]
    FC.RN <- FC.RN[, relevant]
    tp <- sum(abs(FC.TN) == 1 & FC.TN == FC.RN)
    tn <- sum(FC.TN == 0 & FC.TN == FC.RN)
    fp <- sum(FC.TN == 0 & FC.TN != FC.RN)
    fn <- sum(abs(FC.TN) == 1 & FC.TN != FC.RN)
    sens[i, j] <- tp/(tp+fn)
    spec[i, j] <- tn/(tn+fp)
    time[i, j] <- sum(as.numeric(CNOresults[[i]][[j]]$results[, 4]))
  }
}

hyeds.med <- apply(hyeds, 1, median)
sens.med <- apply(sens, 1, median)
spec.med <- apply(spec, 1, median)
time.med <- apply(time, 1, median)
sens.med <- rev(sens.med)
spec.med <- rev(spec.med)
time.med <- rev(time.med)
hyeds.med <- rev(hyeds.med)

lwd <- 2

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_super_nodes.pdf", width = 5, height = 5)

par(mfrow=c(1,1), mar=c(4, 4, 4, 4) + 0.1)

plot(sens.med, ylim = c(0,1), type = "b", xaxt = "n", yaxt = "n", xlab = "nodes", ylab = "sensitivity, specificity", main = "", sub = "", lwd = lwd)
lines(spec.med, type = "b", lty = 2, pch = 2, col = 1, lwd = lwd)
time.med2 <- time.med - min(time.med)
time.med2 <- time.med2/max(time.med2)
#time.med2 <- (time.med2 + 3)/max((time.med2 + 3))
lines(time.med2, type = "b", col = 1, lty = 3, pch = 3, lwd = lwd)
axis(2, seq(0,1,0.1), seq(0,1,0.1))
axis(3, 1:length(hyeds.med), round(hyeds.med), col = "black", col.ticks = "black", col.axis = "black")
axis(1, 1:length(time.med), round(seq(10, 30, (30-10)/4)))
axis(4, seq(0,1,0.1), round(seq(min(time.med), max(time.med), (max(time.med) - min(time.med))/10)/60*12/60, 1), col = "black", col.ticks = "black", col.axis = "black")
mtext("hyper-edges", side=3, line=3, cex.lab=1,las=0, col="black")
mtext("time in hours", side=4, line=3, cex.lab=1,las=0, col="black")
#legend(1, 0.6, legend = c("sensitivity", "specificity", "time"), col = c(1,1,2), lty = 1:3, pch = 1:3, lwd = lwd)

dev.off()

par(ask=T)
for (j in 1:length(CNOinput[[1]])) {
  for (i in 1:length(CNOinput)) {
    print(sum(CNOinput[[i]][[j]]$TN == 1))
    print(length(model$reacID))
    plotDnf(CNOinput[[i]][[j]]$model1$reacID[as.logical(CNOinput[[i]][[j]]$TN)], simulate = list(stimuli = stimuli, inhibitors = NULL))
  }
}

############# other stuff:

source("Boutros10.svn/method/scripts/cnopt.mod.R")

gates <- c("S1+S2=S3", "S1+S2=S4", "S1+S2=S5", "S3+S4+S5=S6", "S4+S5=S7", "!S6=S8") # nem friendly

## gates <- c("S1=S3", "S1=S4", "S1=S5", "S2=S3", "S2=S4", "S2=S5", "S3=S6", "S4=S6", "S5=S6", "S4=S7", "S5=S7", "!S6=S8") # nem unfriendly

## gates <- c("S1+S2=S3", "S1+S2=S4", "S1+S2=S5", "S3=S6", "S4=S6", "S4+S5=S7", "!S6=S8") # combination

## gates <- c("S1=S3", "S2=S4", "S1+S2=S5", "S3=S6", "S4+S5=S7", "!S6+!S7=S8") # incorrect pkn gtn

## gates <- c("S1=S3", "S2=S4", "S2=S5", "S3=S6", "S4=S7", "S5=D3", "S6+S7=S8") # two gtn in one dataset

## gates <- c("S2=S3", "!S2=S3", "S2=S4", "!S2=S4", "S2=S5", "!S2=S5", "S3+S4+S5=S6", "S3+S5=S7", "S6=S8") # two gtn in one dataset

## gates <- c("S1+S2=S3", "S1+S2=S4", "S1+S2=S5", "S3=S6", "S4=S6", "S4+S5=S7", "!S6=S8") # hidden node

TN <- rep(0,length(NCNOcutCompExp$reacID)) # target network

TN[which(model$reacID %in% gates)] <- 1

bString <- TN
NCNOcutCompExp$reacID[which(bString == 1)]

model <- NCNOcutCompExp
##plotBinary(bString, model)
plotDnf(model$reacID[which(bString == 1)], CNOlist = CNOlist)

#redString <- reduceGraph(bString, CNOlist, model) # does more than absorption no
#plotBinary(redString, model)

########################################
######## here data generation for one and highdimensional:   #######
########################################

source("Boutros10.svn/method/scripts/cnopt.mod.R")
bitString <- TN

resECNOlist<-residualError(CNOlist)
Fields4Sim<-prep4sim(NCNOcutCompExp)
SimList <- Fields4Sim
ModelCut <- NCNOcutCompExp
ModelCut$interMat <- ModelCut$interMat[, as.logical(bitString)]
ModelCut$notMat <- ModelCut$notMat[, as.logical(bitString)]
ModelCut$reacID <- ModelCut$reacID[as.logical(bitString)]
SimListCut <- cutSimList(SimList, bitString)
indexList <- indicesNCNOcutComp
SimResults <- simulateStatesRecursive(CNOlist = CNOlist, model = ModelCut, bString = rep(1, length(ModelCut$reacID)))

timePoint <- "t1"
sizeFac <- 1e-04
NAFac <- 1
nInTot <- length(which(NCNOcutCompExp$interMat == -1))
Model <- ModelCut

CNOlist@signals[[2]] <- SimResults

colnames(SimResults) <- colnames(CNOlist@signals[[2]])

rownames(SimResults) <- rownames(CNOlist@signals[[2]])
CNOlist@signals[[2]] <- SimResults

SCompMat <- computeFc(CNOlist, t(CNOlist@signals[[2]]))

####### two gtns:

SimResultsII <- SimResults

SimResults <- cbind(SimResults1, SimResults)

######## weighted SCompMat:

SCompMat <- SCompMat[, colnames(NEMlist$fc)]

par(mfrow=c(1,2))

hist(SCompMat)

for (i in rownames(SCompMat)) {
  SCompMat[i, grep(i, colnames(SCompMat))] <- SCompMat[i, grep(i, colnames(SCompMat))]*2
}

hist(SCompMat)

#### save tns and data:

ANDtn <- TN
ANDdata <- SimResults

ORtn <- TN
ORdata <- SimResults

ORdata2 <- SimResults

## for (i in inhibitors) {
##   SCompMat <- SCompMat[, -grep(paste(i, ".*", i, sep = ""), colnames(SCompMat))]
##   SCompMat <- SCompMat[, -grep(paste("Ctrl.*", i, sep = ""), colnames(SCompMat))]
## }

##heatmapOP(SCompMat, cexCol = 0.6, cexRow = 0.6, Colv = F, Rowv = F, breaks = c(-1,-0.5,0.5,1))

################ pkns with missing edge:

PKN2 <- PKN

### A:
PKN2$reacID <- PKN$reacID[-grep("S4=", PKN$reacID)]
PKN2$notMat <- PKN$notMat[, -grep("S4=", colnames(PKN$notMat))]
PKN2$interMat <- PKN$interMat[, -grep("S4=", colnames(PKN$interMat))]

### B:
PKN2$reacID <- PKN$reacID[-grep("S7=", PKN$reacID)]
PKN2$notMat <- PKN$notMat[, -grep("S7=", colnames(PKN$notMat))]
PKN2$interMat <- PKN$interMat[, -grep("S7=", colnames(PKN$interMat))]

### C:
PKN2$reacID <- PKN$reacID[-grep("S2=S5", PKN$reacID)]
PKN2$notMat <- PKN$notMat[, -grep("S2=S5", colnames(PKN$notMat))]
PKN2$interMat <- PKN$interMat[, -grep("S2=S5", colnames(PKN$interMat))]

checkSignals(CNOlist,PKN2)
indices<-indexFinder(CNOlist,PKN2,verbose=TRUE)
NCNOindices<-findNONC(PKN2,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN2,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate=100)
resECNOlist<-residualError(CNOlist)
Fields4Sim<-prep4sim(NCNOcutCompExp)

modelA <- NCNOcutCompExp

modelB <- NCNOcutCompExp

modelC <- NCNOcutCompExp

modelORG <- model

###### pkn with missing node:

CNOlistORG <- CNOlist
modelORG <- model

sifMatrix <- numeric()
sifMatrix <- rbind(sifMatrix, c("S1", 1, "S3"))
sifMatrix <- rbind(sifMatrix, c("S1", 1, "S5"))
sifMatrix <- rbind(sifMatrix, c("S2", 1, "S3"))
sifMatrix <- rbind(sifMatrix, c("S2", 1, "S5"))
sifMatrix <- rbind(sifMatrix, c("S2", -1, "S3"))
sifMatrix <- rbind(sifMatrix, c("S2", -1, "S5"))
sifMatrix <- rbind(sifMatrix, c("S3", 1, "S6"))
sifMatrix <- rbind(sifMatrix, c("S3", 1, "S7"))
sifMatrix <- rbind(sifMatrix, c("S5", 1, "S6"))
sifMatrix <- rbind(sifMatrix, c("S5", 1, "S7"))
sifMatrix <- rbind(sifMatrix, c("S5", -1, "S6"))
sifMatrix <- rbind(sifMatrix, c("S5", -1, "S7"))
sifMatrix <- rbind(sifMatrix, c("S6", 1, "S8"))
sifMatrix <- rbind(sifMatrix, c("S6", -1, "S8"))
sifMatrix <- rbind(sifMatrix, c("S7", 1, "S8"))
sifMatrix <- rbind(sifMatrix, c("S7", -1, "S8"))

write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN2 <- readSIF("temp.sif")
unlink("temp.sif")

CNOlist <- dummyCNOlist(stimuli, inhibitors[-which(inhibitors %in% "S4")], maxStim = 2, maxInhibit = 2) # length(inhibitors))

checkSignals(CNOlist,PKN2)
indices<-indexFinder(CNOlist,PKN2,verbose=TRUE)
NCNOindices<-findNONC(PKN2,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN2,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate=100)
resECNOlist<-residualError(CNOlist)
Fields4Sim<-prep4sim(NCNOcutCompExp)

model <- NCNOcutCompExp

bString <- absorption(rep(1, length(model$reacID)), model)

plotDnf(model$reacID[as.logical(bString)])

SimResults <- SimResults[-grep("S4", rownames(SimResults)), ]

##### do the runs:

single.run <- T

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source.R")

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source2.R") # two gtns

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source.local.R")

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source.locVga.R")

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source3.R")

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source4.R")

Egenes <- NEMlist$fc
Sgenes <- SCompMat[, colnames(Egenes)]
score <- cor(t(Sgenes), t(Egenes), method = "s")
hist(score)

################## some results stuff for a single or the last run:

tmp <- computeScoreNemT1(CNOlist, model = NCNOcutCompExp, RN, simList = NULL, indexList = NULL, NAFac = 1, NEMlist = NEMlist, tellme = 1, parameters = parameters, sizeFac = sizeFac, method = method)

EtoS <- tmp$EtoS

pdf("sim_hidden.pdf", width = 10, height = 10)

par(mfrow=c(2,2))

plotDnf(modelORG$reacID[-grep("\\+", modelORG$reacID)], CNOlist = CNOlist, legend = 0)
plotDnf(model$reacID[-grep("\\+", model$reacID)], CNOlist = CNOlist, legend = 0)
plotDnf(modelORG$reacID[as.logical(TN)], CNOlist = CNOlist, legend = 0)
plotDnf(model$reacID[as.logical(RN)], CNOlist = CNOlist, legend = 0)

dev.off()

source("Boutros10.svn/method/scripts/cnopt.mod.R")
residuals <- resBNEM(RN, CNOlist, model, NEMlist, parameters, method, cut = T)

resDiff3 <- residuals$resDiff3

pdf("sim_hidden_res2.pdf", width = ncol(resDiff3), height = nrow(resDiff3))

res.breaks <- seq(-max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))), max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))), (max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))) - -max(abs(min(resDiff3, na.rm = T)), abs(max(resDiff3, na.rm = T))))/11)

p1 <- heatmapOP(resDiff3[, 1:(which.max(apply(resDiff3, 2, sum)) - 1)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (positive effects)", sub = "", breaks = res.breaks, colorkey = FALSE, cexCol = 1)

p2 <- heatmapOP(resDiff3[, (which.max(apply(resDiff3, 2, sum)) + 1):ncol(resDiff3)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (negative effects)", sub = "", breaks = res.breaks, colorkey = T, cexCol = 1)

a <- length(1:(which.max(apply(resDiff3, 2, sum)) - 1))
b <- length((which.max(apply(resDiff3, 2, sum)) + 1):ncol(resDiff3))

print(p1, position=c(0, 0, (a/(a+b))-0.1, 1), more=TRUE)
print(p2, position=c((a/(a+b))-0.09, 0, 1, 1))

dev.off()

############### egene heatmaps:

source("Boutros10.svn/method/scripts/cnopt.mod.R")
geneLists <- list()
par(ask=T)
for (i in c(1:ncol(CNOlist@signals[[1]]))) {
  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = RN, Sgene = i, Egenes = 50, parameters=parameters, disc = 0, method = method, dendrogram = "none", cexRow = 0.8, cexCol = 0.4, soft = T, xrot = 90)
}

for (i in 1:length(geneLists)) {
  if (!is.null(dim(geneLists[[i]]))) { next() }
  NEMlist$exprs <- NEMlist$exprs[-which(rownames(NEMlist$exprs) %in% rownames(geneLists[[i]]$genesInfo)), ]
}
  
parameters <- list(cutOffs = c(0.7,0.7,0.25), scoring = c(2,0.5), fitBest = 0) # cutOffs = c(0,0.7, fitBest = 10 seemed to be good
popSize <- 1000 # 100
stallGenMax <- 100
maxGens <- Inf
selection <- c("tournament", 2) # "sus"
initBstring<-rep(0,length(NCNOcutCompExp$reacID)) # teste mal completely disconnected als start netzwerk
initBstring <- bString

Sgenes <- NULL

parallelList <- list()

parallelList[[1]] <- c(0,0,0,2)
parallelList[[2]] <- c("rhskl4", "rhskl5", "rhskl9", "rhskl11")

parallelList <- NULL

source("Boutros10.svn/method/scripts/cnopt.mod.R")
## system.time(
crossRes <- crossTalk(
                      parallel=parallelList,
                      CNOlist=CNOlist,
                      NEMlist=NEMlist,
                      popSize = popSize, # 100
                      maxTime = Inf, # 86400*0.8, # 86400s = 24h
                      maxGens = maxGens, # Inf
                      stallGenMax = stallGenMax, # 10000
                      elitism = ceiling(popSize*0.1), # 10
                      inversion = ceiling(popSize*0.1), # 10
                      verbose=TRUE,
                      selection = selection,
                      parameters = parameters,
                      Sgenes = Sgenes
                      )
## )

dtm.cutoff <- 0.5

dtms <- numeric()

for (i in names(crossRes)) {

  dtms <- c(dtms, crossRes[[i]]$result$dtmRatio)

}

par(mfrow=c(2,3))

for (i in names(crossRes)[order(dtms)]) {
  tmp <- crossRes[[i]]
  print(i)
  print(crossRes[[i]]$result$dtmRatio)
  if (crossRes[[i]]$result$dtmRatio >= dtm.cutoff) { next() }
  if (max(tmp$result$bString) == 0 | sum(tmp$result$bString[3:4]) == 2) {
    plot(1)
  } else {
    plotBinary(tmp$result$bString, tmp$model)
    
  }
}

##############r

source("Kube011_BCR_CD40_Inhibitoren/publication/3.data.simulation.source.R")

bString <- RN

source("Boutros10.svn/method/scripts/cnopt.mod.R")
computeScoreNemT1(CNOlist, model = model, bString, simList = NULL, indexList = NULL, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

plotDnf(model$reacID[as.logical(bString)])

geneLists <- list()
source("Boutros10.svn/method/scripts/cnopt.mod.R")
par(ask=T)
for (i in c(1:(ncol(CNOlist@signals[[1]])-0))) {
  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = bString, Sgene = i, Egenes = 50, parameters=parameters, disc = 0, method = method, dendrogram = "none", cexRow = 0.8, cexCol = 0.7, soft = T, sub = "")
}
