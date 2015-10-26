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
negative <- 1

sifMatrix <- numeric()

for (i in stimuli) {
  for (j in inhibitors[1:6]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1 & negative) {
      direction <- "-1"
    } else {
      direction <- "1"
    }
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(j, direction, i))
    } else {
      sifMatrix <- rbind(sifMatrix, c(i, direction, j))
    }
  }
}
for (i in inhibitors[1:6]) {
  for (j in inhibitors[7:12]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1 & negative) {
      direction <- "-1"
    } else {
      direction <- "1"
    }
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(j, direction, i))
    } else {
      sifMatrix <- rbind(sifMatrix, c(i, direction, j))
    }
  }
}
for (j in inhibitors[7:12]) {
  for (k in inhibitors[13:18]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1 & negative) {
      direction <- "-1"
    } else {
      direction <- "1"
    }
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(k, direction, j))
    } else {
      sifMatrix <- rbind(sifMatrix, c(j, direction, k))
    }
  }
}
for (k in inhibitors[13:18]) {
  for (l in inhibitors[19:24]) {
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1 & negative) {
      direction <- "-1"
    } else {
      direction <- "1"
    }
    if (sample(c(0,1), 1, prob=c(0.9, 0.1)) == 1) {
      sifMatrix <- rbind(sifMatrix, c(l, direction, k))
    } else {
      sifMatrix <- rbind(sifMatrix, c(k, direction, l))
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

hierarchy <- list()

hierarchy[[1]] <- paste("S", 1:6, sep = "")

for (i in 1:4) {
  hierarchy[[(i+1)]] <- paste("I", i, 1:6, sep = "")
}

if (negative) {
  pdf("dissertation/diss/Thesis/gfx/BCR_super_pkn_neg_cyclic.pdf", height = 10, width = 13)
}else {
  pdf("dissertation/diss/Thesis/gfx/BCR_super_pkn_cyclic.pdf", height = 10, width = 13)
}

plotDnf(model$reacID[-grep("\\+", model$reacID)], legend = 0, width = 1.25, edgecol = edgecol, stimuli = stimuli, height = 1)#, hierarchy = hierarchy)

dev.off()

length(model$reacID)

if (negative) {
  cluster.file <- "cluster/nodesCneg.RData"
} else {
  cluster.file <- "cluster/nodesC.RData"
}

save.image(cluster.file)

####################### athene:

hpc.nodes <- 1
hpc.ppn <- 8
walltime <- "48:00:00"
pmem <- "1000"

noises <- c(0.1,0.25,0.5,0.5,1,2)
types <- c(rep("disc", 3), rep("cont", 3))

for (i in 1:length(noises)) {

  write(paste("#!/bin/bash\n#PBS -l walltime=", walltime, "\n#PBS -l pmem=", pmem, "mb\n#PBS -l nodes=", hpc.nodes, ":ppn=", hpc.ppn, "\ncd /spang.compdiag/user/pim28150/\nRBioCscript cluster/nodesC.R 'data=", cluster.file, "' 'noise=", noises[i], "' 'type=", types[i], "'", sep = ""), file = paste("cluster/nodesC", i, ".sh", sep = ""))

}

dev.off()

system("ssh -X athene", intern = F)

## qsub /spang.compdiag/user/pim28150/cluster/nodesC1.sh -l walltime=00:3:00 -q express -o /spang.compdiag/user/pim28150/cluster/temp.out -e /spang.compdiag/user/pim28150/cluster/temp.err

## qsub /spang.compdiag/user/pim28150/cluster/temp.sh -q spang -o /spang.compdiag/user/pim28150/cluster/temp.out -e /spang.compdiag/user/pim28150/cluster/temp.err

for i in `seq 1 6`; do

  qsub /spang.compdiag/user/pim28150/cluster/nodesC$i.sh -q serial -o /spang.compdiag/user/pim28150/cluster/nodesC$i.out -e /spang.compdiag/user/pim28150/cluster/nodesC$i.err

  rm -f /spang.compdiag/user/pim28150/cluster/nodesC$i.sh
  sleep 5
  done

exit

print("start of the script")

library(methods)
source("github/trunk/method/cnopt.mod.R")
library(CellNOptR)

args <- commandArgs()
file <- gsub("data=", "", args[grep("data=", args)])
noise <- as.numeric(gsub("noise=", "", args[grep("noise=", args)]))
type <- gsub("type=", "", args[grep("type=", args)])

load(file)

print("everything loaded")

######################## do fundamental super pkn sim:

source("Boutros10.svn/method/scripts/cnopt.mod.R")

method <- "s"
parameters <- list(cutOffs = c(0,1,0), scoring = c(0.1,0.9,2))
sizeFac <- 10^-10

verbose <- TRUE
popSize <- 100
stallGenMax <- 10
maxGens <- Inf

parallel <- 8#list(4, "rhskl4")

count <- 0

CNOinput <- list()

CNOresults <- list()

TN.list <- list()

for (n in rev(round(seq(10, 30, (30-10)/4)))) {
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
      if (length(nodes.tmp[-which(nodes.tmp %in% core)]) < n - length(core)) {
        core <- c(core, nodes.tmp[-which(nodes.tmp %in% core)])
      } else {
        core <- c(core, sample(nodes.tmp[-which(nodes.tmp %in% core)], n - length(core)))
      }
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
    if ("disc" %in% type) {
      noisy <- sample(1:length(NEMlist$exprs), floor(noise*length(NEMlist$exprs)))
      NEMlist$exprs[noisy] <- 1 - NEMlist$exprs[noisy]
    }
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
    if ("cont" %in% type) {
      NEMlist$fc <- NEMlist$fc + rnorm(length(NEMlist$fc), 0, noise)
    }
    start <- Sys.time()
    CNOresults[[count]][[c]] <- gaBinaryNemT1(parallel=parallel,CNOlist=CNOlist,NEMlist=NEMlist,model=model1,initBstring=initBstring,popSize = popSize, maxTime = Inf, maxGens = maxGens, stallGenMax = stallGenMax, elitism = ceiling(popSize*0.01),inversion = ceiling(popSize*0.01),verbose=TRUE,parameters = parameters,sizeFac = sizeFac,method = method)
    end <- Sys.time()
    print(end - start)
  }
}

## save(CNOlist, CNOresults, CNOinput, file = paste("superpkn_nodes_cycles_", noise, "_", paste(type, collapse = "_"), ".RData", sep = ""))

source("temp.R")

load(paste("superpkn_nodes_cycles_", noise, "_", paste(type, collapse = "_"), ".RData", sep = ""))

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

pdf(paste("dissertation/diss/Thesis/gfx/super_nodes_cyclic_", type, "_", noise, ".pdf", sep = ""), width = 5, height = 5)

par(mfrow=c(1,1), mar=c(4, 4, 4, 4) + 0.1)

plot(sens.med, ylim = c(0,1), type = "b", xaxt = "n", yaxt = "n", xlab = "nodes", ylab = "sensitivity, specificity", main = "", sub = "", lwd = lwd)
lines(spec.med, type = "b", lty = 2, pch = 2, col = 1, lwd = lwd)
time.med2 <- time.med - min(time.med)
time.med2 <- time.med2/max(time.med2)
#time.med2 <- (time.med2 + 3)/max((time.med2 + 3))
#lines(time.med2, type = "b", col = 1, lty = 3, pch = 3, lwd = lwd)
axis(2, seq(0,1,0.1), seq(0,1,0.1))
axis(3, 1:length(hyeds.med), round(hyeds.med), col = "black", col.ticks = "black", col.axis = "black")
axis(1, 1:length(time.med), round(seq(10, 30, (30-10)/4)))
#axis(4, seq(0,1,0.1), round(seq(min(time.med), max(time.med), (max(time.med) - min(time.med))/10)/60*12/60, 1), col = "black", col.ticks = "black", col.axis = "black")
mtext("hyper-edges", side=3, line=3, cex.lab=1,las=0, col="black")
#mtext("time in hours", side=4, line=3, cex.lab=1,las=0, col="black")
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
