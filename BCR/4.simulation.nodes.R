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

## PKN <- readSIF("Kube011_BCR_CD40_Inhibitoren/publication/PKNjulio.sif")
## stimuli <- c("il6", "tgfa", "ins", "tnf", "il1")
## inhibitors <- unique(gsub("!", "", unlist(strsplit(PKN$reacID, "="))))
## inhibitors <- inhibitors[-which(inhibitors %in% stimuli)]

plotDnf(PKN$reacID)
length(PKN$reacID)

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

plotDnf(model$reacID[-grep("\\+", model$reacID)], legend = 0, width = 0.5, edgecol = edgecol)

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

err <- list()

TN.list <- list()

for (n in rev(round(seq(10, 30, (30-10)/4)))) {
  ## n <- 30
  count <- count + 1
  CNOinput[[count]] <- list()
  CNOresults[[count]] <- list()
  err[[count]] <- numeric()
  for (c in 1:10) {
    ## c <- 1
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
    if (n == 30) {
      TN <- sample(c(0,1), length(model1$reacID), replace = T)
      if (length(grep(paste(paste("S", 1:6, sep = ""), collapse = "|"), model1$reacID[as.logical(TN)])) == 0) {
        TN[sample(grep(paste(stimuli[which(stimuli %in% nodes)], collapse = "|"), model1$reacID), sample(1:length(grep(paste(stimuli[which(stimuli %in% nodes)], collapse = "|"), model1$reacID)), 1))] <- 1
      }
      TN <- reduceGraph(TN, model1, CNOlist)
      TN.list[[c]] <- model1$reacID[as.logical(TN)]
    } else {
      TN <- numeric(length(model1$reacID))
      TN[which(model1$reacID %in% TN.list[[c]])] <- 1
    }
    SimResults <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = TN)
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
    start <- Sys.time()
    CNOresults[[count]][[c]] <- gaBinaryNemT1(parallel=parallel,CNOlist=CNOlist,NEMlist=NEMlist,model=model1,initBstring=initBstring,popSize = popSize, maxTime = Inf, maxGens = maxGens, stallGenMax = stallGenMax, elitism = ceiling(popSize*0.01),inversion = ceiling(popSize*0.01),verbose=TRUE,parameters = parameters,sizeFac = sizeFac,method = method)
    end <- Sys.time()
    print(end - start)
    RN <- CNOresults[[count]][[c]]$bString
    ERS.RN <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = RN)
    ERS.TN <- simulateStatesRecursive(CNOlist = CNOlist, model = model1, bString = TN)
    err[[count]][c] <- sum(ERS.RN - ERS.TN != 0)
    print(paste("error: ", err[[count]][c], " normalized error: ", err[[count]][c]/length(ERS.TN), sep = ""))
  }
}

source("temp.R")

computeScoreNemT1(CNOlist, model = model, TN, simList = NULL, indexList = NULL, NAFac = 1, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

par(mfrow=c(1,2))
plotDnf(model1$reacID[as.logical(TN)], CNOlist = CNOlist, legend = 0)
plotDnf(model1$reacID[as.logical(RN)], CNOlist = CNOlist, legend = 0)

## save(CNOlist, CNOresults, CNOinput, file = "superpkn_nodes_new.RData")

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
time.med2 <- (time.med2 + 3)/max((time.med2 + 3))
lines(time.med2, type = "b", col = 2, lty = 3, pch = 3, lwd = lwd)
axis(2, seq(0,1,0.1), seq(0,1,0.1))
axis(3, 1:length(hyeds.med), round(hyeds.med), col = "red", col.ticks = "red", col.axis = "red")
axis(1, 1:length(time.med), round(seq(10, 30, (30-10)/4)))
axis(4, seq(0.75,1,0.025), round(seq(min(time.med), max(time.med), (max(time.med) - min(time.med))/10)/60*12/60, 1), col = "red", col.ticks = "red", col.axis = "red")
mtext("hyper-edges", side=3, line=3, cex.lab=1,las=0, col="red")
mtext("time in hours", side=4, line=3, cex.lab=1,las=0, col="red")
#legend(1, 0.6, legend = c("sensitivity", "specificity", "time"), col = c(1,1,2), lty = 1:3, pch = 1:3, lwd = lwd)

dev.off()
