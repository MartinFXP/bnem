############################

X11.options(type="Xlib")

library(limma)
source("Boutros10.svn/method/scripts/cnopt.mod.R")
library(CellNOptR)
library(annotate)
library(hgu133plus2.db)

load("Kube011_BCR_CD40_Inhibitoren/publication/NEMlist.Combat.RData")

#################### check it out:

colnames(NEMlist$exprs)[-grep("CD40|Vivit|Myc|JSH", colnames(NEMlist$exprs))] # 72 = BCR*6 + Ctrl*6 + BCR*(6 + 4) + Ctrl*(6 + 4)

################### start:

batches <- paste("Batch", 1:6, sep = "")

runs <- paste("Exp", 1:2, sep = "")

stimuli <-  c("BCR", "CD40")
stimuli <- stimuli[1]
stimuli.not <- c("BCR", "CD40")
stimuli.not <- stimuli.not[-which(stimuli.not %in% stimuli)]

inhibitors <- c("Erk1", "Erk2", "Ikk2", "Jnk", "p38", "Myc", "Pi3k", "Tak1")

inhibitors <- inhibitors[-which(inhibitors %in% "Erk2")]
inhibitors <- inhibitors[-which(inhibitors %in% "Myc")]
inhibitors <- gsub("Erk1", "Erk", inhibitors)

if (stimuli %in% "CD40") {
  inhibitors <- inhibitors[-which(inhibitors %in% c("Erk","Pi3k"))]
}

## inhibitors <- inhibitors[-which(inhibitors %in% c("Erk", "Pi3k", "Tak1"))]
signals <- inhibitors

signals <- c(stimuli, inhibitors)

hierarchy <- list(c("BCR"), c("PI3K"), c("TAK1"), c("IKK2", "P38", "JNK"))

## signals <- c(stimuli, inhibitors, paste("M", 1:20, sep = ""))

control <- c("Neg", "Ctrl")
times <- "none"

colnames(NEMlist$exprs) <- gsub("Erk1", "Erk", gsub("xErk2", "", colnames(NEMlist$exprs)))
colnames(NEMlist$fc) <- gsub("Erk1", "Erk", gsub("_Erk2", "", colnames(NEMlist$fc)))

CNOlist <- cnoFromData(NEMlist$exprs, times = times, stimuli = stimuli, inhibitors = inhibitors, signals = signals)

## columns einschränken:

for (i in inhibitors) {
  if (length(grep(paste(i, ".*", i, sep = ""), colnames(NEMlist$fc))) > 0) {
    NEMlist$fc <- NEMlist$fc[, -grep(paste(i, ".*", i, sep = ""), colnames(NEMlist$fc))]
  }
}

NEMlist$fc <- NEMlist$fc[, -grep("Myc", colnames(NEMlist$fc))]

dim(NEMlist$fc)

## take out not used stimuli:

if (length(stimuli.not) > 0) {

  NEMlist$fc <- NEMlist$fc[, -grep(paste(stimuli.not, collapse = "|"), colnames(NEMlist$fc))]

}

dim(NEMlist$fc)

############### get vivit out of there (what is vivit anyway?):

NEMlist$fc <- NEMlist$fc[, -grep("Vivit", colnames(NEMlist$fc))]

dim(NEMlist$fc)

############# remove AFFYX controls:

affy.ind <- grep("AFFX", rownames(NEMlist$exprs))

NEMlist$exprs <- NEMlist$exprs[-affy.ind, ]
NEMlist$fc <- NEMlist$fc[-affy.ind, ]
NEMlist$pvals <- NEMlist$pvals[-affy.ind, ]
NEMlist$B <- NEMlist$B[-affy.ind, ]

dim(NEMlist$fc)

########## B-NEM:

## reconstruction prior:

## sifMatrix <- numeric()
## for (i in inhibitors) {
##   if (i %in% c("p38", "Ikk2", "Jnk")) {
##     sifMatrix <- rbind(sifMatrix, c(i, 1, paste("M", 1:9, sep = "")))
##   }
## }
## for (i in inhibitors) {
##   for (j in stimuli) {
##     sifMatrix <- rbind(sifMatrix, c(j, 1, i, rep("", ncol(sifMatrix) - 3)))
##   }
##   for (j in inhibitors) {
##     if (i %in% j) { next() }
##     sifMatrix <- rbind(sifMatrix, c(j, 1, i, rep("", ncol(sifMatrix) - 3)))
##   }
## }

## hierarchical prior:

sifMatrix <- numeric()
for (i in stimuli) {
    sifMatrix <- rbind(sifMatrix, c(i, 1, inhibitors))
}
sifMatrix <- rbind(sifMatrix, c("Pi3k", 1, inhibitors[-which(inhibitors %in% "Pi3k")], rep("", ncol(sifMatrix) - 2 - length(inhibitors[-which(inhibitors %in% "Pi3k")]))))
sifMatrix <- rbind(sifMatrix, c("Tak1", 1, inhibitors[-which(inhibitors %in% c("Pi3k", "Tak1"))], rep("", ncol(sifMatrix) - 2 - length(inhibitors[-which(inhibitors %in% c("Pi3k", "Tak1"))]))))

## moduls:

if (length(grep("M[1-9]|M[1-2][0-9]", signals)) > 0) {
  for (i in inhibitors[-which(inhibitors %in% c("Pi3k", "Tak1"))]) {
    for (j in signals[-which(signals %in% c(stimuli, inhibitors))]) {
      sifMatrix <- rbind(sifMatrix, c(i, 1, j, rep("", ncol(sifMatrix) - 3)))
      ##sifMatrix <- rbind(sifMatrix, c(i, -1, j, rep("", ncol(sifMatrix) - 3)))
    }
  }
}

write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

stimuli <- stimuli[which(stimuli %in% PKN$namesSpecies)]
inhibitors <- inhibitors[which(inhibitors %in% PKN$namesSpecies)]
signals <- signals[which(signals %in% PKN$namesSpecies)]

checkSignals(CNOlist,PKN)
indices<-indexFinder(CNOlist,PKN,verbose=TRUE)
NCNOindices<-findNONC(PKN,indices,verbose=TRUE)
NCNOcut<-cutNONC(PKN,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate = 100)
resECNOlist<-residualError(CNOlist)
Fields4Sim<-prep4sim(NCNOcutCompExp)

model <- NCNOcutCompExp

bString <- numeric(length(model$reacID))
bString[-grep("\\+", model$reacID)] <- 1
bString[grep(paste(paste("M", 1:100, sep = ""), collapse = "|"), model$reacID)] <- 0

source("Boutros10.svn/method/scripts/cnopt.mod.R")

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_pkn.pdf", height = 5, width = 5)

graph <- model$reacID[-grep("\\+", model$reacID)]

graph <- gsub("Tak1", "TAK1", gsub("Erk", "ERK", gsub("p38", "P38", gsub("Ikk2", "IKK2", gsub("Jnk", "JNK", gsub("Pi3k", "PI3K", graph))))))

if (length(grep("M[0-9]|M[0-9][0-9]", graph)) > 0) {
  plotDnf(graph[-grep("M[0-9]|M[0-9][0-9]", graph)], stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", length(graph[-grep("M[0-9]|M[0-9][0-9]", graph)])))
} else {
  plotDnf(graph, stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", length(graph)))
}

dev.off()

length(model$reacID)

################# backup data:

source("Boutros10.svn/method/scripts/cnopt.mod.R")

## NEMlist.backup <- NEMlist

############### pca to get egenes?:

NEMlist <- NEMlist.backup

pca.model <- prcomp(t(NEMlist$fc))

pdf("temp.pdf", height = 10, width = 10)

biplot(pca.model, col = c("black", "transparent"))

dev.off()

library(rgl)
plotPCA <- function(x, nGroup) {
    n <- ncol(x) 
    if(!(n %in% c(2,3))) { # check if 2d or 3d
        stop("x must have either 2 or 3 columns")
    }

    fit <- hclust(dist(x), method="complete") # cluster
    groups <- cutree(fit, k=nGroup)

    if(n == 3) { # 3d plot
        plot3d(x, col=groups, type="s", size=1, axes=F)
        axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
        grid3d("x")
        grid3d("y")
        grid3d("z")
    } else { # 2d plot
        maxes <- apply(abs(x), 2, max)
        rangeX <- c(-maxes[1], maxes[1])
        rangeY <- c(-maxes[2], maxes[2])
        plot(x, col=groups, pch=19, xlab=colnames(x)[1], ylab=colnames(x)[2], xlim=rangeX, ylim=rangeY)
        lines(c(0,0), rangeX*2)
        lines(rangeY*2, c(0,0))
    }
}

plotPCA(pca.model$x[,1:2], 2)
plotPCA(pca.model$x[,1:3], 2)

NEMlist$fc <- (t(pca.model$x%*%t(pca.model$rotation)) + pca.model$center) # "original"

NEMlist$fc <- t(pca.model$x)

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", inhibitors, sep = ""), collapse = "|"), colnames(NEMlist$fc))]

NEMlist$fc <- NEMlist$fc[1:9, ]

dim(NEMlist$fc)

############# reduce egenes:

NEMlist <- NEMlist.backup

NEMlist <- dataRed(NEMlist, stimuli, c(1), receptors = NULL, inhibitors, c(log2(1.5),0,0,log2(1.5)), direction = NULL)

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", inhibitors, sep = ""), collapse = "|"), colnames(NEMlist$fc))]

dim(NEMlist$fc)

################ get top scoring genes genes:

NEMlist <- NEMlist.backup

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", inhibitors, sep = ""), collapse = "|"), colnames(NEMlist$fc))]

## NEMlist$fc <- NEMlist$fc[, -grep("Ikk2.*p38|Ikk2.*Jnk|Jnk.*p38", colnames(NEMlist$fc))]

dim(NEMlist$fc)

corsel <- 0.01

############# get top correlated genes with pkn:

pkn.targets <- NULL

source("Boutros10.svn/method/scripts/cnopt.mod.R")

bString <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)
bString[grep(paste(paste("M", 1:999, sep = ""), collapse = "|"), model$reacID)] <- 0

method <- "s"
parameters <- list(cutOffs = c(0,1,0), scoring = c(0.1,0.9,2))

tmp <- computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 1, parameters = parameters, method = method)

EtoS <- tmp$EtoS

pkn.targets <- which(rownames(NEMlist$fc) %in% rownames(EtoS)[which(EtoS[, 4] <= quantile(EtoS[, 4], corsel))])

############ get top correlated genes with nem:

nem.targets <- NULL

bString <- rep(1, length(model$reacID))
bString[grep(paste(paste("M", 1:999, sep = ""), collapse = "|"), model$reacID)] <- 0
bString <- absorptionII(bString, model)

method <- "s"
parameters <- list(cutOffs = c(0,1,0), scoring = c(0.1,0.9,2))

tmp <- computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 1, parameters = parameters, method = method)

EtoS <- tmp$EtoS

nem.targets <- which(rownames(NEMlist$fc) %in% rownames(EtoS)[which(EtoS[, 4] <= quantile(EtoS[, 4], corsel))])

########### combine:

NEMlist <- NEMlist.backup

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", inhibitors, sep = ""), collapse = "|"), colnames(NEMlist$fc))]

## NEMlist$fc <- NEMlist$fc[, -grep("Ikk2.*p38|Ikk2.*Jnk|Jnk.*p38", colnames(NEMlist$fc))]

dim(NEMlist$fc)

NEMlist$fc <- NEMlist$fc[unique(c(pkn.targets, nem.targets)), ]

NEMlist$exprs <- NEMlist$exprs[unique(c(pkn.targets, nem.targets)), ]

dim(NEMlist$fc)

########### intersect:

NEMlist <- NEMlist.backup

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", inhibitors, sep = ""), collapse = "|"), colnames(NEMlist$fc))]

## NEMlist$fc <- NEMlist$fc[, -grep("Ikk2.*p38|Ikk2.*Jnk|Jnk.*p38", colnames(NEMlist$fc))]

dim(NEMlist$fc)

NEMlist$fc <- NEMlist$fc[intersect(pkn.targets, nem.targets), ]

NEMlist$exprs <- NEMlist$exprs[intersect(pkn.targets, nem.targets), ]

dim(NEMlist$fc)

####################### cnopipeline:

method <- "s"
parameters <- list(cutOffs = c(0,1,0), scoring = c(0.1,0.9,2))
sizeFac <- 10^-10

verbose <- TRUE
popSize <- 100
stallGenMax <- 10
maxGens <- Inf
initBstring <- reduceGraph(rep(0,length(model$reacID)), model)

parallel <- 12#list(8, "rhskl11")

#### prep for cluster:

## NEMlist1 <- NEMlist

## cluster.file <- paste("cluster/B-NEM/BCR", 1, ".RData", sep = "")

## save(NEMlist1, CNOlist, model, parameters, method, sizeFac, initBstring, parallel, verbose, popSize, stallGenMax, file = cluster.file)

## write(paste("#!/bin/bash\n#PBS -l walltime=01:00:00\n#PBS -l pmem=500mb\n#PBS -l nodes=", 1, ":ppn=", 8, "\ncd /spang.compdiag/user/pim28150/\nRBioCscript cluster/B-NEM/gaDyn.R 'run=", 1, "' 'data=", cluster.file, "'", sep = ""), file = paste("cluster/B-NEM/BCRga", 1, ".sh", sep = ""))

## dev.off()

## system("ssh -X athene", intern = F)

## qsub /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.sh -l walltime=00:01:00 -q express -o /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.out -e /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.err

## qsub /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.sh -q serial -o /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.out -e /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.err
## rm -f /spang.compdiag/user/pim28150/cluster/B-NEM/BCRga1.sh

##### not on cluster:

for (i in 1:10) {

source("Boutros10.svn/method/scripts/cnopt.mod.R")
CNOresult <- gaBinaryNemT1(
               parallel=parallel,
               CNOlist=CNOlist,
               NEMlist=NEMlist,
               model=model,
               initBstring=initBstring,
               popSize = popSize, # 100
               maxTime = Inf, # 86400*0.8, # 86400s = 24h
               maxGens = maxGens, # Inf
               stallGenMax = stallGenMax, # 10000
               elitism = ceiling(popSize*0.01),
               inversion = ceiling(popSize*0.01),
               verbose=TRUE,
               parameters = parameters,
               sizeFac = sizeFac,
               method = method
               )

#### get from cluster:

save(CNOresult, NEMlist, CNOlist, model, parameters, sizeFac, popSize, file = paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", paste(stimuli, collapse = "_"), "_moduls_result", format(Sys.time(), "%Y.%m.%d.%H:%m:%S"), ".RData", sep = ""))

}

CNOruns <- list()

scores <- numeric()

files <- character()

count <- 0

egenes <- 2629 # 2456 for the moduls (egenes based on priors) and 750 for the core bcr 1 and 2629 for core bcr log2(1.5)

modellen <- 32

for (i in list.files("Kube011_BCR_CD40_Inhibitoren/publication/results/")) {
  if (length(grep(paste(paste(stimuli, collapse = "_"), "_moduls_result", sep = ""), i)) > 0) {
    popSize <- 1
    load(paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", i, sep = ""))
    if (popSize != 100) {
      unlink(i)
      next()
    }
    if (dim(NEMlist$fc)[1] != egenes | length(model$reacID) != modellen) {
      next()
    }
    count <- count + 1
    CNOruns[[count]] <- CNOresult
    files[count] <- paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", i, sep = "")
    scores[count] <- computeScoreNemT1(CNOlist, model = model, CNOresult$bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)
    print(scores[count])
    print(dim(NEMlist$fc))
  }
}
load(files[which.min(scores)])

## Rprof(NULL)
## summaryRprof(filename="Rprof.out") # simulatorT1 takes the most time (not anymore?)

drawScores(CNOresult)
#dev.print(paste(paste(stimuli, collapse = "_"), "_moduls_search.pdf", sep = ""), device = pdf)

#colSums(CNOresult$StringsTol)/nrow(CNOresult$StringsTol)
initSeed <- bString <- CNOresult$bString
scaffold.size <- nrow(CNOresult$stringsTol)
scaffold <- apply(CNOresult$stringsTol, 2, sum)
scaffold <- scaffold/scaffold.size
comp <- cbind(NCNOcutCompExp$reacID,scaffold, CNOresult$bString)
comp[which(as.numeric(comp[, 2]) > 0.5 |  as.numeric(comp[, 3]) == 1), ]
CNOresultTmp2 <- CNOresult

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_res_freq.pdf", height = 5, width = 5)

best.lty <- rep("dashed", sum(scaffold > 0.5 | bString == 1))

best.lty[which(model$reacID[which(scaffold > 0.5 | bString == 1)] %in% model$reacID[as.logical(bString)])] <- "solid"

source("Boutros10.svn/method/scripts/cnopt.mod.R")

graph <- model$reacID[which(scaffold > 0.5 | bString == 1)]

graph <- gsub("Tak1", "TAK1", gsub("Erk", "ERK", gsub("p38", "P38", gsub("Ikk2", "IKK2", gsub("Jnk", "JNK", gsub("Pi3k", "PI3K", graph))))))

pkn.graph <- c("BCR=PI3K", "PI3K=TAK1", "TAK1=IKK2", "TAK1=P38", "TAK1=JNK", "TAK1=ERK")

if (sum(pkn.graph %in% graph) > 0) {
  pkn.graph <- pkn.graph[-which(pkn.graph %in% graph)]
}

g <- plotDnf(c(pkn.graph, graph), CNOlist = CNOlist, freq = c(rep(0, length(pkn.graph)), scaffold[which(scaffold > 0.5 | bString == 1)]), width = 1, nodecol = "white", bordercol = "black", lty = c(rep("solid", length(pkn.graph)), best.lty, "dashed", "dashed"), labelcol = "transparent", edgelwd = 2, edgecol = c(rep("transparent", length(pkn.graph)), rep("black", length(graph)+2)))

dev.off()

system.time(
computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)
)

############ get the rel tol networks from all runs:

scaffold.size <- 0
scaffold <- numeric(ncol(CNOresult$stringsTol))
giant.stringsTol <- numeric(ncol(CNOresult$stringsTol))

for (i in 1:length(CNOruns)) {
  scaffold.size <- scaffold.size + nrow(CNOruns[[i]]$stringsTol) 
  scaffold <- scaffold + apply(CNOruns[[i]]$stringsTol, 2, sum)
  giant.stringsTol <- rbind(giant.stringsTol, CNOruns[[i]]$stringsTol)
}

colnames(giant.stringsTol) <- model$reacID

scaffold <- scaffold/scaffold.size

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_res_freq_full.pdf", height = 5, width = 5)

freq.cut <- 0.5

best.lty <- rep("dashed", sum(scaffold > freq.cut | bString == 1))

best.lty[which(model$reacID[which(scaffold > freq.cut | bString == 1)] %in% model$reacID[as.logical(bString)])] <- "solid"

source("Boutros10.svn/method/scripts/cnopt.mod.R")

g <- plotDnf(model$reacID[which(scaffold > freq.cut | bString == 1)], CNOlist = CNOlist, freq = scaffold[which(scaffold > freq.cut | bString == 1)], edgecol = "black", width = 1, nodecol = "white", bordercol = "black", lty = best.lty, edgelwd = 2.5, labelcol = "transparent")

dev.off()

mut.ex <- matrix(0, length(model$reacID[which(scaffold > 0.5 | bString == 1)]), length(model$reacID[which(scaffold > 0.5 | bString == 1)]))
rownames(mut.ex) <- model$reacID[which(scaffold > 0.5 | bString == 1)]
colnames(mut.ex) <- model$reacID[which(scaffold > 0.5 | bString == 1)]
               
for (i in model$reacID[which(scaffold > 0.5 | bString == 1)]) {
  for (j in model$reacID[which(scaffold > 0.5 | bString == 1)]) {
    mut.ex[which(rownames(mut.ex) %in% i), which(rownames(mut.ex) %in% j)] <- sum(giant.stringsTol[, which(colnames(giant.stringsTol) %in% i)] == 1 & giant.stringsTol[, which(colnames(giant.stringsTol) %in% j)] == 1)/(sum(giant.stringsTol[, which(colnames(giant.stringsTol) %in% i)] == 1) + sum(giant.stringsTol[, which(colnames(giant.stringsTol) %in% j)] == 1) - sum(giant.stringsTol[, which(colnames(giant.stringsTol) %in% i)] == 1 & giant.stringsTol[, which(colnames(giant.stringsTol) %in% j)] == 1))
  }
}

rownames(mut.ex) <- gsub("=", ' -> ', rownames(mut.ex))

rownames(mut.ex)[length(rownames(mut.ex))] <- "(BCR, Pi3k) -> Jnk"

colnames(mut.ex) <- rownames(mut.ex)

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_mutex.pdf", height = 7, width = 8)

heatmapOP(mut.ex, Colv = F, Rowv = F, breaks = seq(0, 1, 0.01), main = "", sub = "", xrot = 45)

dev.off()

############### predictive accuracy (learn zeta):

source("Boutros10.svn/method/scripts/cnopt.mod.R")

CNOresults <- list()

targets <- list()

count <- 0

zetas <- 10^(0:-10)

parallel <- 12#list(8, "rhskl11")

for (i in 1:10) {

  NEMlist1 <- NEMlist

  targets[[i]] <- sample(1:nrow(NEMlist1$fc), floor(0.5*nrow(NEMlist1$fc)))

  NEMlist1$fc <- NEMlist1$fc[targets[[i]], ]
  NEMlist1$exprs <- NEMlist1$exprs[targets[[i]], ]

  dim(NEMlist1$fc)

  for (i in  zetas) {

    print(paste("zeta: ", i, sep = ""))

    count <- count + 1
    CNOresults[[count]] <- gaBinaryNemT1(
                             parallel=parallel,
                             CNOlist=CNOlist,
                             NEMlist=NEMlist1,
                             model=model,
                             initBstring=initBstring,
                             popSize = popSize, # 100
                             maxTime = Inf, # 86400*0.8, # 86400s = 24h
                             maxGens = maxGens, # Inf
                             stallGenMax = stallGenMax, # 10000
                             elitism = ceiling(popSize*0.01),
                             inversion = ceiling(popSize*0.01),
                             verbose=TRUE,
                             parameters = parameters,
                             sizeFac = i,
                             method = method
                             )

  }

}

## save(CNOresults, targets, file = "Kube011_BCR_CD40_Inhibitoren/publication/BCR_zeta_train.RData")

load("Kube011_BCR_CD40_Inhibitoren/publication/BCR_zeta_train.RData")

pred <- matrix(0, length(zetas), 10)

mods <- matrix(0, length(zetas), 10)

graph.size <- matrix(0, length(zetas), 10)

node.num <- matrix(0, length(zetas), 10)

count <- 0

for (i in 1:10) {

  NEMlist1 <- NEMlist

  NEMlist1$fc <- NEMlist1$fc[targets[[i]], ]
  NEMlist1$exprs <- NEMlist1$exprs[targets[[i]], ]

  dim(NEMlist1$fc)

  NEMlist2 <- NEMlist

  NEMlist2$fc <- NEMlist2$fc[-targets[[i]], ]
  NEMlist2$exprs <- NEMlist2$exprs[-targets[[i]], ]

  dim(NEMlist2$fc)

  for (j in 1:length(zetas)) {

    count <- count + 1
    pred[j, i] <- 1 - computeScoreNemT1(CNOlist, model = model, CNOresults[[count]]$bString, NEMlist = NEMlist2, parameters = parameters, sizeFac = zetas[j]*0, method = method)
    graph.tmp <- model$reacID[as.logical(CNOresults[[count]]$bString)]
    nodes.tmp <- unique(unlist(strsplit(unlist(strsplit(graph.tmp, "=")), "\\+")))
    mods[j, i] <- length(grep("M[0-9]|M[0-9][0-9]", nodes.tmp))
    graph.size[j, i] <- length(unlist(strsplit(graph.tmp, "\\+")))
    node.num[j, i] <- length(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(graph.tmp, "=")), "\\+")))))

  }

}

pred.mat <- pred
graph.size.mat <- graph.size
node.num.mat <- node.num

heatmapOP(pred.mat, Colv = F, Rowv = F, breaks = seq(min(pred.mat), max(pred.mat), 0.001))
heatmapOP(graph.size.mat, Colv = F, Rowv = F, breaks = seq(min(graph.size.mat), max(graph.size.mat), 0.5))
heatmapOP(node.num.mat, Colv = F, Rowv = F, breaks = seq(min(node.num.mat), max(node.num.mat), 0.5))

stat.fun <- median

pred <- apply(pred.mat, 1, stat.fun)
graph.size <- apply(graph.size.mat, 1, stat.fun)
node.num <- apply(node.num.mat, 1, stat.fun)

lwd <- 2

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_zeta.pdf", width = 5, height = 5)

par(mfrow=c(1,1), mar=c(4, 4, 4, 4) + 0.1)

pred2 <- pred - min(pred)
mods2 <- mods - min(mods)
mods2 <- (mods2/max(mods2))*max(pred2)
graph.size2 <- graph.size - min(graph.size)
node.num2 <- node.num - min(graph.size)
node.num2 <- (node.num2/max(graph.size2))*max(pred2)
graph.size2 <- (graph.size2/max(graph.size2))*max(pred2)

## this without moduls:

plot(pred2, type = "b", xaxt = "n", yaxt = "n", main = "", ylab = "score of independent dataset", xlab = "zeta values", lwd = lwd)
lines(graph.size2, type = "b", lty = 2, pch = 2, col = "red", lwd = lwd)
lines(node.num2, type = "b", lty = 3, pch = 3, col = "blue", lwd = lwd)
axis(2, seq(min(pred2), max(pred2), (max(pred2) - min(pred2))/9), round(seq(min(pred), max(pred), (max(pred) - min(pred))/9), 3))
axis(1, 1:11, zetas)
axis(4, seq(min(graph.size2), max(graph.size2), (max(graph.size2) - min(graph.size2))/(max(graph.size) - min(graph.size))), seq(min(graph.size), max(graph.size), 1), col = "red", col.ticks = "red", col.axis = "red")
mtext("graph size / number of nodes", side=4, line=3, cex.lab=1, las=0, col="red")

dev.off()

## this with moduls:

## plot(pred2, type = "b", xaxt = "n", yaxt = "n", main = "", ylab = "score of independent dataset", xlab = "zeta values")
## lines(mods2, type = "b", lty = 2, pch = 2, col = "red")
## axis(2, seq(min(pred2), max(pred2), (max(pred2) - min(pred2))/9), round(seq(min(pred), max(pred), (max(pred) - min(pred))/9), 3))
## axis(1, 1:11, zetas)
## axis(4, seq(min(mods2), max(mods2), (max(mods2) - min(mods2))/(max(mods) - min(mods))), seq(min(mods), max(mods), 1), col = "red", col.ticks = "red", col.axis = "red")
## mtext("modules used", side=4, line=3, cex.lab=1, las=0, col="red")

## dev.off()

## save(CNOresults, targets, zetas, NEMlist, file = "BCR_zeta_train.RData")

############################################ local search

initSeed <- reduceGraph(rep(0, length(model$reacID)), model, CNOlist)

initSeed[grep("=M", model$reacID)] <- 0

initSeed <- bString

par(ask=F)
source("Boutros10.svn/method/scripts/cnopt.mod.R")
system.time(
  localString <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist, model=model, parameters=parameters, verbose = TRUE, parallel=parallel, parallel2=1, seeds=1, initSeed = initSeed, sizeFac = sizeFac, method = method, draw = TRUE, max.steps = Inf)
  )

## drawLocal(localString)

## graph <- model$reacID[as.logical(localString$bStrings[1, ])]

## plotDnf(graph, stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", length(unlist(strsplit(graph, "\\+")))+1))

## ## NEMlist$fc[which(mget(rownames(NEMlist$fc), hgu133plus2SYMBOL) %in% "AICDA"), ]

## computeScoreNemT1(CNOlist, model = NCNOcutCompExp, localString$bStrings[1, ], NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

## bString <- localString$bStrings[1, ]

##### check on original e-genes:

computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

geneLists <- list()
source("Boutros10.svn/method/scripts/cnopt.mod.R")
par(ask=F)
for (i in 1:(ncol(CNOlist@signals[[1]]))) {

  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = bString, Sgene = i, Egenes = 30, method = method, parameters=parameters, cexRow = 0.8, cexCol = 1, EtoS = TRUE, affyIds = T, soft = T, sub = "", sizeFac = sizeFac, disc = 0, Colv = T, Rowv = T, dendrogram = "col", csc = FALSE, xrot = 60, aspect = "iso")

}

for (i in 1:length(geneLists)) {

  pdf(paste("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_", colnames(CNOlist@signals[[1]])[i], ".pdf", sep = ""), height = 15, width = 10)

  data <- geneLists[[i]]$data

  rownames(data) <- gsub("Tak1", "TAK1", gsub("Erk", "ERK", gsub("p38", "P38", gsub("Ikk2", "IKK2", gsub("Jnk", "JNK", gsub("Pi3k", "PI3K", rownames(data)))))))
  
  colnames(data) <- paste(c(rep("(BCR+)", 10), "(Control)"), paste("(", gsub(",vs,", ") vs (", gsub("BCR", "BCR+", gsub("Erk", "ERK-", gsub("Ikk2", "IKK2-", gsub("Jnk", "JNK-", gsub("p38", "P38-", gsub("_", ",", gsub("Tak1", "TAK1-", gsub("Pi3k", "PI3K-", colnames(data)))))))))), ")", sep = ""), sep = " vs ")

  print(heatmapOP(data, cexRow = 0.8, cexCol = 1, Colv = T, Rowv = F, dendrogram = "col", csc = FALSE, xrot = 45, aspect = "iso", breaks = seq(-1, 1, 0.05), main = "", sub = ""))

  dev.off()
  
}

i <- 3
pdf(paste("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_", colnames(CNOlist@signals[[1]])[i], "_grey.pdf", sep = ""), height = 8, width = 7)

geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = bString, Sgene = i, Egenes = 25, method = method, parameters=parameters, cexRow = 0.8, cexCol = 1, EtoS = TRUE, affyIds = T, soft = T, sub = "", sizeFac = sizeFac, disc = 0, Colv = T, Rowv = T, dendrogram = "col", csc = FALSE, xro = 60, col = "Greys")

dev.off()

i <- 13
pdf(paste(Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/"BCR_", colnames(CNOlist@signals[[1]])[i], "_grey.pdf", sep = ""), height = 8, width = 7)

geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = bString, Sgene = i, Egenes = 25, method = method, parameters=parameters, cexRow = 0.8, cexCol = 1, EtoS = TRUE, affyIds = T, soft = T, sub = "", sizeFac = sizeFac, disc = 0, Colv = T, Rowv = T, dendrogram = "col", csc = FALSE, xro = 60, col = "Greys")

dev.off()

################## checking on egenes starts here:

computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_res.pdf", height = 5, width = 5)

graph <- model$reacID[which(bString == 1)]

plotDnf(graph, stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", length(unlist(strsplit(graph, "\\+")))+1))

## plotDnf(graph[-grep("M[0-9]|M[0-9][0-9]", graph)], stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", length(unlist(strsplit(graph[-grep("M[0-9]|M[0-9][0-9]", graph)], "\\+")))+1))

dev.off()

EtoS <- computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 1, parameters = parameters, sizeFac = sizeFac, method = method)

moduldist <- numeric()

for (i in 1:ncol(CNOlist@signals[[1]])) {
  moduldist <- c(moduldist, sum(EtoS$EtoS[, 2] == i))
}

modulgraph <- model$reacID[which(bString == 1)]

for (i in grep("M[0-9]|M[0-9[0-9]", colnames(CNOlist@signals[[1]]))) {
  if (sum(EtoS$EtoS[, 2] == i) > quantile(moduldist, 0.5)) {
    modulgraph <- c(modulgraph, paste(colnames(CNOlist@signals[[1]])[i], "=", sum(EtoS$EtoS[, 2] == i), sep = ""))
  } else {
    if (length(grep(colnames(CNOlist@signals[[1]])[i], modulgraph)) > 0) {
      modulgraph <- modulgraph[-grep(paste("=", colnames(CNOlist@signals[[1]])[i], "$|", colnames(CNOlist@signals[[1]])[i], "=", sep = ""), modulgraph)]
    }
  }
}

pdf("BCR_res_moduls.pdf", height = 10, width = 15)

plotDnf(modulgraph, legend = 0, width = 1, edgecol = "black")

dev.off()

## check residuals:

source("Boutros10.svn/method/scripts/cnopt.mod.R")
residuals <- resBNEM(bString, CNOlist, model, NEMlist, parameters, method, sizeFac = sizeFac, cut = F)

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_residuals.pdf", width = 10, height = 5)

max(abs(min(residuals$resDiff3, na.rm = T)), abs(max(residuals$resDiff3, na.rm = T)))

res.breaks <- seq(-max(abs(min(residuals$resDiff3, na.rm = T)), abs(max(residuals$resDiff3, na.rm = T))), max(abs(min(residuals$resDiff3, na.rm = T)), abs(max(residuals$resDiff3, na.rm = T))), (max(abs(min(residuals$resDiff3, na.rm = T)), abs(max(residuals$resDiff3, na.rm = T))) - -max(abs(min(residuals$resDiff3, na.rm = T)), abs(max(residuals$resDiff3, na.rm = T))))/100)

res.tmp <- residuals$resDiff3
colnames(res.tmp) <- gsub("BCR_vs_", "", colnames(res.tmp))

p1 <- heatmapOP(res.tmp[, 1:((ncol(res.tmp)-1)/2)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (positive effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = FALSE, aspect = "iso")

p2 <- heatmapOP(res.tmp[, (((ncol(res.tmp)-1)/2) + 2):ncol(res.tmp)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (negative effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = T, aspect = "iso")

print(p1, position=c(0, 0, .45, 1), more=TRUE)
print(p2, position=c(.45, 0, 1, 1))

dev.off()

#### add hyperedge and relearn:

CNOlist2 <- addSignal("M1", CNOlist)

model2 <- addEdge(c("Erk=M", "Ikk2=M1", "p38=M1", "Jnk=M1", "Tak1=M1", "Pi3k=M1"), CNOlist2, model)

#### add manual hyperedge (no relearn):

model2 <- addEdge(c("Erk=Pi3k", "Erk=Ikk2"), CNOlist, model)#addEdge(c("Ikk2=Erk", "p38=Erk", "Jnk=Erk", "Tak1=Erk"), CNOlist, model)

bString2 <- numeric(length(model2$reacID))
        
model2$reacID

bString2[which(model2$reacID %in% c("BCR=Erk", "Erk=Pi3k", "Erk=Ikk2", "BCR=Tak1", "Pi3k=Jnk", "Pi3k=p38", "Tak1+Jnk=Erk", "Tak1+Ikk2+p38=Erk"))] <- 1

#bString2[which(model2$reacID %in% c("BCR=Ikk2", "BCR=Pi3k", "BCR=Tak1", "Pi3k=Jnk", "Pi3k=p38", "Tak1+Jnk=Erk", "Tak1+Ikk2+p38=Erk"))] <- 1

plotDnf(model2$reacID[which(bString2 == 1)], CNOlist = CNOlist)

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_res_debug.pdf", height = 3, width = 3)

graph2 <- model2$reacID[which(bString2 == 1)]

#graph2 <- graph2[-grep("M9|M[1-7]|M[0-9][0-9]|Erk|Ikk2|Jnk", graph2)]

#plotDnf(graph2, legend = 0, width = 1, edgecol = c(rep("black", 5), "red", rep("black", 2), "red", rep("black", 1)))

plotDnf(graph2, legend = 0, width = 1, edgecol = "black")

dev.off()

residualsII <- resBNEM(bString2, CNOlist, model2, NEMlist, parameters, method)

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_residualsII.pdf", width = 10, height = 5)

max(abs(min(residualsII$resDiff3, na.rm = T)), abs(max(residualsII$resDiff3, na.rm = T)))

res.breaks <- seq(-max(abs(min(residualsII$resDiff3, na.rm = T)), abs(max(residualsII$resDiff3, na.rm = T))), max(abs(min(residualsII$resDiff3, na.rm = T)), abs(max(residualsII$resDiff3, na.rm = T))), (max(abs(min(residualsII$resDiff3, na.rm = T)), abs(max(residualsII$resDiff3, na.rm = T))) - -max(abs(min(residualsII$resDiff3, na.rm = T)), abs(max(residualsII$resDiff3, na.rm = T))))/100)

res.tmp <- residualsII$resDiff3
colnames(res.tmp) <- gsub("BCR_vs_", "", colnames(res.tmp))

p1 <- heatmapOP(res.tmp[, 1:((ncol(res.tmp)-1)/2)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (positive effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = FALSE, aspect = "iso")

p2 <- heatmapOP(res.tmp[, (((ncol(res.tmp)-1)/2) + 2):ncol(res.tmp)], bordercol = "grey", Colv = F, Rowv = F, main = "residuals (negative effects)", sub = "", xrot = "60", breaks = res.breaks, colorkey = T, aspect = "iso")

print(p1, position=c(0, 0, .45, 1), more=TRUE)
print(p2, position=c(.45, 0, 1, 1))

dev.off()

computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method) - computeScoreNemT1(CNOlist, model = model2, bString2, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method) # positive says score improved

pdf("temp.pdf", width = 10, height = 5)

########### relearn network:

for (i in 1:10) {

initBstring <- numeric(length(model2$reacID))
initBstring[which(model2$reacID %in% model$reacID[as.logical(bString)])] <- 1

source("Boutros10.svn/method/scripts/cnopt.mod.R")
CNOresult2 <- gaBinaryNemT1(
               parallel=parallel,
               CNOlist=CNOlist2,
               NEMlist=NEMlist,
               model=model2,
               initBstring=initBstring,
               popSize = popSize, # 100
               maxTime = Inf, # 86400*0.8, # 86400s = 24h
               maxGens = maxGens, # Inf
               stallGenMax = stallGenMax, # 10000
               elitism = ceiling(popSize*0.01),
               inversion = ceiling(popSize*0.01),
               verbose=TRUE,
               parameters = parameters,
               sizeFac = sizeFac,
               method = method
               )

save(CNOresult2, NEMlist, CNOlist2, model2, parameters, sizeFac, popSize, file = paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", paste(stimuli, collapse = "_"), "_moduls_result2", format(Sys.time(), "%Y.%m.%d.%H:%m:%S"), ".RData", sep = ""))

}

scores <- numeric()

files <- character()

count <- 0

egenes <- 2629 # 2456 # 2456 for the moduls (egenes based on priors) and 750 for the core bcr 1 and 2629 for core bcr log2(1.5)

for (i in list.files("Kube011_BCR_CD40_Inhibitoren/publication/results/")) {
  if (length(grep(paste(paste(stimuli, collapse = "_"), "_moduls_result2", sep = ""), i)) > 0) {
    popSize <- 1
    load(paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", i, sep = ""))
    if (popSize != 100) {
      unlink(i)
      next()
    }
    if (dim(NEMlist$fc)[1] != egenes) {
      next()
    }
    count <- count + 1
    files[count] <- paste("Kube011_BCR_CD40_Inhibitoren/publication/results/", i, sep = "")
    scores[count] <- computeScoreNemT1(CNOlist2, model = model2, CNOresult2$bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)
    print(scores[count])
    print(dim(NEMlist$fc))
  }
}
load(files[which.min(scores)])

initSeed <- bString2 <- CNOresult2$bString

## local:

initSeed <- reduceGraph(rep(1, length(model2$reacID)), model2, CNOlist2)

initSeed[grep("=M", model2$reacID)] <- 0

par(ask=F)
source("Boutros10.svn/method/scripts/cnopt.mod.R")
system.time(
  localString2 <- localSearch(CNOlist=CNOlist2, NEMlist=NEMlist, model=model2, parameters=parameters, verbose = TRUE, parallel=parallel, parallel2=1, seeds=1, initSeed = initSeed, sizeFac = sizeFac, method = method, draw = TRUE, max.steps = Inf)
  )

bString2 <- localString2$bStrings[1, ]

## check improvements with modul:

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_res_debug.pdf", height = 8, width = 5)

graph2 <- model2$reacID[which(bString2 == 1)]

#graph2 <- graph2[-grep("M9|M[1-7]|M[0-9][0-9]|Erk|Ikk2|Jnk", graph2)]

#plotDnf(graph2, legend = 0, width = 1, edgecol = c(rep("black", 5), "red", rep("black", 2), "red", rep("black", 1)))

plotDnf(graph2, legend = 0, width = 1, edgecol = "black")

dev.off()

computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method) - computeScoreNemT1(CNOlist2, model = model2, bString2, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method) # positive says score improved

residualsIII <- resBNEM(bString2, CNOlist2, model2, NEMlist, parameters, method)

##### check on new e-gene distribution:

computeScoreNemT1(CNOlist2, model = model2, bString2, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

geneLists <- list()
source("Boutros10.svn/method/scripts/cnopt.mod.R")
par(ask=T)
for (i in 1:(ncol(CNOlist2@signals[[1]]))) {

  pdf(paste("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_", colnames(CNOlist2@signals[[1]])[i], "_debug.pdf", sep = ""), height = 8, width = 7)

  geneLists[[i]] <- validateGraph(CNOlist2, NEMlist, model = model2, bString = bString2, Sgene = i, Egenes = 50, method = method, parameters=parameters, cexRow = 0.8, cexCol = 1, EtoS = TRUE, affyIds = T, soft = T, sub = "", sizeFac = sizeFac, disc = 0, Colv = T, Rowv = T, dendrogram = "col", csc = FALSE, xrot = 60)

  dev.off()
  
}

EtoS2 <- computeScoreNemT1(CNOlist2, model = model2, bString2, NEMlist = NEMlist, tellme = 1, parameters = parameters, sizeFac = sizeFac, method = method)

#### fisher.test for enrichment:

genes <- list(M1 = unlist(mget(rownames(EtoS2$EtoS)[which(EtoS2$EtoS[, 2] == which(colnames(CNOlist2@signals[[1]]) %in% "M1"))], hgu133plus2SYMBOL)), Global = unlist(mget(rownames(NEMlist$fc), hgu133plus2SYMBOL)))

load("mgsa_broad_genesymbols/GSEAlistII.RData")

parallel.gsea <- 8

gsea.set <- GSEAlistII[[5]]#c(GSEAlistII[[2]], GSEAlistII[[3]], GSEAlistII[[5]])

gsea.res <- myGsea(genes, gsea.set, parallel = parallel.gsea, conservative = T)

gsea.res2 <- gsea.res[which(gsea.res$test %in% "M1"), ]

gsea.res2[order(gsea.res2$qvals)[1:30], ]

##### plot each modul connection for new residual network:

for (i in 1:20) {

  if (sum(EtoS2$EtoS[, 2] == which(colnames(CNOlist@signals[[1]]) %in% paste("M", i, sep = ""))) == 0) {
    bString2[grep(paste("M", i, sep = ""), model2$reacID)] <- 0
  }
  
  graph.tmp <- model2$reacID[which(bString2 == 1)]

  modul.tmp <- graph.tmp[grep(paste("M", i, "$", sep = ""), graph.tmp)]

  graph.tmp <- graph.tmp[-grep("M[1-9]|M[1-9][0-9]", graph.tmp)]

  graph.tmp <- c(graph.tmp, modul.tmp)

  graph.tmp <- c(graph.tmp, paste("M", i, "=", sum(EtoS2$EtoS[, 2] == which(colnames(CNOlist@signals[[1]]) %in% paste("M", i, sep = ""))), sep = ""))

  pdf(paste("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_M", i, "_graph.pdf", sep = ""), height = 8, width = 5)

  edgecol.tmp <- rep("black", length(unlist(strsplit(graph.tmp, "\\+")))+length(unlist(strsplit(graph.tmp, "\\+")))-length(graph.tmp))

  if (length(grep("!", unlist(strsplit(graph.tmp, "\\+")))) > 0) {
    edgecol.tmp[grep("!", unlist(strsplit(graph.tmp, "\\+")))+(2:length(grep("\\+", graph.tmp))-1)] <- "red"
  }

  plotDnf(graph.tmp, legend = 0, width = 1, edgecol = edgecol.tmp, stimuli = "BCR")

  dev.off()

}

### do equivalent graph:

pdf("Kube011_BCR_CD40_Inhibitoren/publication/pdf/bioinformatics/gfx/BCR_resII.pdf", height = 5, width = 5)

bStringE <- bString
bStringE[which(model$reacID %in% "Pi3k=Tak1")] <- 0
bStringE[which(model$reacID %in% "Tak1=p38")] <- 0
bStringE[which(model$reacID %in% "Tak1=Erk")] <- 0
bStringE[which(model$reacID %in% "BCR=Tak1")] <- 1
bStringE[which(model$reacID %in% "Pi3k+Tak1=p38")] <- 1
bStringE[which(model$reacID %in% "Pi3k+Tak1=Erk")] <- 1

graph <- model$reacID[which(bStringE == 1)]

plotDnf(graph[-grep("M[0-9]|M[0-9][0-9]", graph)], stimuli = "BCR", legend = 0, width = 1, edgecol = rep("black", 10))

dev.off()

#### end

