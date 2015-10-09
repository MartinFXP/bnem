############################

X11.options(type="Xlib")

library(limma)
source("github/trunk/method/cnopt.mod.R")
library(CellNOptR)
library(annotate)
library(hgu133plus2.db)

load("publications/BNEM/NEMlist.Combat.RData")

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

signals <- c(stimuli, inhibitors)

hierarchy <- list(c("BCR"), c("PI3K"), c("TAK1"), c("IKK2", "P38", "JNK"))

control <- c("Neg", "Ctrl")
times <- "none"

colnames(NEMlist$exprs) <- gsub("Erk1", "Erk", gsub("xErk2", "", colnames(NEMlist$exprs)))
colnames(NEMlist$fc) <- gsub("Erk1", "Erk", gsub("_Erk2", "", colnames(NEMlist$fc)))

CNOlist <- cnoFromData(NEMlist$exprs, times = times, stimuli = stimuli, inhibitors = inhibitors, signals = signals)

## get columns:

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

############### get vivit out of there:

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

## recon prior for ikk2, p38, jnk and no feed forward:

negation <- F
sifMatrix <- numeric()
for (i in stimuli) {
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

source("github/trunk/method/cnopt.mod.R")

pdf("publications/BNEM/BCR_pkn.pdf", height = 5, width = 5)

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

source("github/trunk/method/cnopt.mod.R")

## NEMlist.backup <- NEMlist

############# reduce egenes:

source("github/trunk/method/cnopt.mod.R")

NEMlist <- NEMlist.backup

NEMlist <- dataRed(NEMlist, stimuli, rep(1, length(stimuli)), receptors = NULL, inhibitors, c(log2(2),0,log2(1.5),0), direction = -1)

NEMlist$fc <- NEMlist$fc[, -grep(paste(paste("Ctrl_vs_", c(inhibitors, "Pi3k", "Erk"), sep = ""), collapse = "|"), colnames(NEMlist$fc))]

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

parallel <- list(8, "rhskl11")

############### predictive accuracy (learn zeta):

source("github/trunk/method/cnopt.mod.R")

CNOresults <- list()

targets <- list()

count <- 0

zetas <- c(1, rev(seq(10^-5,1,0.2)^2), 0) # zetas <- c(10^((-0):(-10), 0) # zetas <- rev(seq(0,1,0.2)^10)

parallel <- list(8, "rhskl11")

trainruns <- 100

for (i in 1:trainruns) {

  NEMlist1 <- NEMlist

  targets[[i]] <- sample(1:nrow(NEMlist1$fc), floor(0.5*nrow(NEMlist1$fc)))

  NEMlist1$fc <- NEMlist1$fc[targets[[i]], ]
  NEMlist1$exprs <- NEMlist1$exprs[targets[[i]], ]

  dim(NEMlist1$fc)

  for (j in zetas) {

    print(paste("run: ", i, sep = ""))

    print(paste("zeta: ", j, sep = ""))

    ## source("github/trunk/method/cnopt.mod.R"); bString <- absorptionII(sample(c(1,0), 10, replace = T), model); computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist1, parameters = parameters, sizeFac = j, method = method)

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
                             sizeFac = j,
                             method = method
                             )

  }

}

## save(CNOresults, targets, zetas, file = "publications/BNEM/BCR_zeta_train_recon.RData")

load("publications/BNEM/BCR_zeta_train_recon.RData")

source("github/trunk/method/cnopt.mod.R")

pred <- matrix(0, length(zetas), trainruns)

mods <- matrix(0, length(zetas), trainruns)

graph.size <- matrix(0, length(zetas), trainruns)

node.num <- matrix(0, length(zetas), trainruns)

count <- 0

for (i in 1:trainruns) {

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
    pred[j, i] <- computeScoreNemT1(CNOlist, model = model, CNOresults[[count]]$bString, NEMlist = NEMlist2, parameters = parameters, sizeFac = zetas[j]*0, method = method, opt = "max")
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

stat.fun <- mean

pred <- apply(pred.mat, 1, stat.fun)
graph.size <- apply(graph.size.mat, 1, stat.fun)
node.num <- apply(node.num.mat, 1, stat.fun)

graph.size <- (graph.size/length(model$reacID))*100
node.num <- (node.num/ncol(CNOlist@signals[[1]]))*100

lwd <- 2

pdf("publications/BNEM/BCR_zeta.pdf", width = 6, height = 6)

par(mfrow=c(1,1), mar=c(4, 4, 4, 4) + 0.1)

pred2 <- pred - min(pred)
mods2 <- mods - min(mods)
mods2 <- (mods2/max(mods2))*max(pred2)
node.num2 <- (node.num/100)*max(pred2)
graph.size2 <- (graph.size/100)*max(pred2)

log <- "y"
if (nchar(log) > 0) {
  pred2 <- pred2 + 1
  node.num2 <- node.num2 + 1
  graph.size2 <- graph.size2 + 1
}

## this without moduls:

plot(pred2, type = "b", xaxt = "n", yaxt = "n", main = "", ylab = "score of test dataset", xlab = "zeta values", lwd = lwd, log = log)
lines(graph.size2, type = "b", lty = 2, pch = 2, col = "black", lwd = lwd)
#lines(node.num2, type = "b", lty = 3, pch = 3, col = "black", lwd = lwd)
axis(2, seq(min(pred2), max(pred2), (max(pred2) - min(pred2))/9), round(seq(min(pred), max(pred), (max(pred) - min(pred))/9), 4))
zetas2 <- c("1", "0.64", "0.36", "0.16", "0.04", "10e-10", "0")
axis(1, 1:length(zetas), zetas2)
##axis(1, seq(0, length(zetas), 2), c(0, zetas[(1:6)*2]))
axis(4, seq(min(pred2), max(pred2), (max(pred2) - min(pred2))/9), round(seq(0,100, ((100 - 0)/(length(seq(min(pred2), max(pred2), (max(pred2) - min(pred2))/9)) - 1))), -1), col = "black", col.ticks = "black", col.axis = "black")
mtext("graph size in %", side=4, line=3, cex.lab=1, las=0, col="black")
#mtext("graph size, number of connected S-genes in %", side=4, line=3, cex.lab=1, las=0, col="black")

dev.off()

##### not on cluster:

for (i in 1:10) {

  source("github/trunk/method/cnopt.mod.R")
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
                 elitism = ceiling(popSize*0.1),
                 inversion = ceiling(popSize*0.1),
                 verbose=TRUE,
                 parameters = parameters,
                 sizeFac = sizeFac,
                 method = method
                 )

  save(CNOresult, NEMlist, CNOlist, model, parameters, sizeFac, popSize, file = paste("publications/BNEM/results/", paste(stimuli, collapse = "_"), "_moduls_result", format(Sys.time(), "%Y.%m.%d.%H:%m:%S"), ".RData", sep = ""))

}

CNOruns <- list()

scores <- numeric()

files <- character()

count <- 0

egenes <- 602

modellen <- 54

for (i in list.files("publications/BNEM/results/")) {
  if (length(grep(paste(paste(stimuli, collapse = "_"), "_moduls_result", sep = ""), i)) > 0) {
    popSize <- 1
    load(paste("publications/BNEM/results/", i, sep = ""))
    if (dim(NEMlist$fc)[1] != egenes | length(model$reacID) != modellen) {
      next()
    }
    count <- count + 1
    CNOruns[[count]] <- CNOresult
    files[count] <- paste("publications/BNEM/results/", i, sep = "")
    scores[count] <- computeScoreNemT1(CNOlist, model = model, CNOresult$bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)
    print(scores[count])
    print(dim(NEMlist$fc))
  }
}
load(files[which.min(scores)])

relTol <- 0.1
initSeed <- bString <- CNOresult$bString
scaffold.size <- nrow(CNOresult$stringsTol[which(CNOresult$stringsTolScores <= (min.score + relTol*min.score)), ])
scaffold <- apply(CNOresult$stringsTol[which(CNOresult$stringsTolScores <= (min.score + relTol*min.score)), ], 2, sum)
scaffold <- scaffold/scaffold.size
comp <- cbind(model$reacID,scaffold, CNOresult$bString)
comp[which(as.numeric(comp[, 2]) > 0.5 |  as.numeric(comp[, 3]) == 1), ]
CNOresultTmp2 <- CNOresult

min.score <- computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

plotDnf(model$reacID[as.logical(bString)])

pdf("publications/BNEM/BCR_res.pdf", height = 7, width = 7)

#plotDnf(c("BCR=PI3K", "BCR=TAK1", "PI3K=OR1", "TAK1=OR1", "OR1=IKK2", "TAK1=ERK", "PI3K=JNK", "PI3K=AND", "IKK2=AND", "JNK=OR2", "AND=OR2", "OR2=P38"), stimuli = "BCR", nodewidth = list(BCR=1,PI3K=1,TAK1=1,IKK2=1,JNK=1,P38=1,ERK=1,OR1=0.5,OR2=0.5,AND=0.5), nodeheight = list(BCR=1,PI3K=1,TAK1=1,IKK2=1,JNK=1,P38=1,ERK=1,OR1=0.5,OR2=0.5,AND=0.5), edgecol = "black", nodecol = list(AND = "grey", OR1 = "grey", OR2="grey"), bordercol = list(AND = "grey", OR1 = "grey", OR2="grey"), nodelabel = list(OR1 = "OR", OR2 = "OR"))

plotDnf(c("BCR=PI3K", "BCR=TAK1", "PI3K=OR1", "TAK1=OR1", "OR1=IKK2", "TAK1=ERK", "PI3K=JNK", "PI3K=AND", "IKK2=AND", "JNK=OR2", "AND=OR2", "OR2=P38", "PI3K=Transcription", "TAK1=Transcription", "IKK2=Transcription", "JNK=Transcription", "P38=Transcription", "ERK=Transcription"), stimuli = "BCR", nodewidth = list(BCR=1,PI3K=1,TAK1=1,IKK2=1,JNK=1,P38=1,ERK=1,OR1=0.5,OR2=0.5,AND=0.5,Transcription=5), nodeheight = list(BCR=1,PI3K=1,TAK1=1,IKK2=1,JNK=1,P38=1,ERK=1,OR1=0.5,OR2=0.5,AND=0.5,Transcription=1), edgecol = c(rep("black", 12), rep("grey", 6)), nodecol = list(AND = "grey", OR1 = "grey", OR2="grey",Transcription="grey"), bordercol = list(AND = "grey", OR1 = "grey", OR2="grey",Transcription="grey"), nodelabel = list(OR1 = "OR", OR2 = "OR"), edgestyle = c(rep("solid", 12), rep("dashed", 6)))

dev.off()

initSeed <- bString <- CNOresult$bString
computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

## initSeed <- bString <- CNOresult$bString
## bString[grep("^Jnk=p38|^Pi3k\\+Ikk2=p38", model$reacID)] <- 0
## bString[grep("^Ikk2\\+Jnk=p38", model$reacID)] <- 1

## computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

##### check on e-genes:

computeScoreNemT1(CNOlist, model = model, bString, NEMlist = NEMlist, tellme = 0, parameters = parameters, sizeFac = sizeFac, method = method)

geneLists <- list()
source("github/trunk/method/cnopt.mod.R")
par(ask=F)
for (i in 1:(ncol(CNOlist@signals[[1]]))) {

  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = bString, Sgene = i, Egenes = 30, method = method, parameters=parameters, cexRow = 0.8, cexCol = 1, EtoS = TRUE, affyIds = F, soft = T, sub = "", sizeFac = sizeFac, disc = 0, Colv = T, Rowv = T, dendrogram = "col", csc = FALSE, xrot = 60, aspect = "iso", breaks = seq(-2,2,0.1))
  #dev.print("temp.pdf", device = pdf)

}

library(biomaRt)

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

mart.att <- listAttributes(mart)

mart.filt <- listFilters(mart)

genes.meta <- getBM(filters= "affy_hg_u133_plus_2", 
                    attributes= c("ensembl_gene_id", "hgnc_symbol", "affy_hg_u133_plus_2"),
                    values= rownames(NEMlist$fc),
                    mart= mart)

par(ask=F)
for (i in 1:length(geneLists)) {

  pdf(paste("publications/BNEM/BCR_", colnames(CNOlist@signals[[1]])[i], ".pdf", sep = ""), height = 15, width = 10)

  data <- geneLists[[i]]$data

  rownames(data) <- gsub("Tak1", "TAK1", gsub("Erk", "ERK", gsub("p38", "P38", gsub("Ikk2", "IKK2", gsub("Jnk", "JNK", gsub("Pi3k", "PI3K", rownames(data)))))))
  
  colnames(data) <- paste(c(rep("(BCR+)", 10), "(Control)"), paste("(", gsub(",vs,", ") vs (", gsub("BCR", "BCR+", gsub("Erk", "ERK-", gsub("Ikk2", "IKK2-", gsub("Jnk", "JNK-", gsub("p38", "P38-", gsub("_", ",", gsub("Tak1", "TAK1-", gsub("Pi3k", "PI3K-", colnames(data)))))))))), ")", sep = ""), sep = " vs ")

  gene.symbols <- character((nrow(data)-2))
  
  for (j in 1:(nrow(data)-2)) {
    gene.tmp <- genes.meta[which(genes.meta[, 3] %in% rownames(data)[j]), 2]
    if (sum(gene.tmp %in% "") > 0) {
      gene.tmp <- gene.tmp[-which(gene.tmp %in% "")]
    }
    gene.tmp <- paste(gene.tmp, collapse = ", ")
    if (nchar(gene.tmp) == 0) {
      gene.tmp <- "NA"
    }
    gene.symbols[j] <- gene.tmp
    ##genes.meta[which(genes.meta[, 3] %in% rownames(data)[1:(nrow(data)-2)]), 2]
  }
  rownames(data)[1:(nrow(data)-2)] <- mget(rownames(data)[1:(nrow(data)-2)], hgu133plus2SYMBOL)
  rownames(data)[which(rownames(data) %in% "NA")] <- gene.symbols[which(rownames(data) %in% "NA")]
  #rownames(data)[1:(nrow(data)-2)] <- gene.symbols
  
  print(heatmapOP(data, cexRow = 0.8, cexCol = 1, Colv = T, Rowv = F, dendrogram = "col", csc = FALSE, xrot = 45, aspect = "iso", main = "", sub = "", breaks = seq(-2, 2, 0.1)))

  dev.off()
  
}

####################### NEM comparison:

source("github/trunk/method/cnopt.mod.R")
library(nem)

## do cutoff learning:

nem.data <- abs(NEMlist$fc[, grep(paste(paste(".*BCR_", inhibitors, "$", sep = ""), collapse = "|"), colnames(NEMlist$fc))])

colnames(nem.data) <- gsub("BCR_vs_BCR_", "", colnames(nem.data))

nem.data <- disc(nem.data, log2(1.5))

nem.res <- nem(nem.data)

pdf("temp.pdf", width = 5, height = 100)
heatmapOP(nem.data, aspect = "iso")
dev.off()

plot(nem.res)

tmp <- graph2adj(nem.res$graph)

tmp <- transitive.reduction(tmp)

tmp1 <- adj2dnf(tmp)

count <- 0

nodelabel <- list()
nodecol <- list()
bordercol <- list()
nodeheight <- list()
nodewidth <- list()

for (i in gsub(".*=", "=", tmp1)) {
  if (length(grep(i, tmp1)) > 1) {
    count <- count + 1
    tmp1 <- gsub(i, paste("=AND", count, sep = ""), tmp1)
    tmp1 <- c(tmp1, paste("AND", count, i, sep = ""))
    nodelabel <- c(nodelabel, "")
    nodecol <- c(nodecol, "transparent")
    bordercol <- c(bordercol, "transparent")
    nodeheight <- c(nodeheight, 0)
    nodewidth <- c(nodewidth, 0)
  }
}

names(nodelabel) <- paste("AND", 1:count, sep = "")
names(nodecol) <- paste("AND", 1:count, sep = "")
names(bordercol) <- paste("AND", 1:count, sep = "")
names(nodeheight) <- paste("AND", 1:count, sep = "")
names(nodewidth) <- paste("AND", 1:count, sep = "")

pdf("publications/BNEM/BCR_nem_res.pdf", height = 5, width = 3)

plotDnf(tmp1, width = 1, nodelabel = nodelabel, nodeheight = nodeheight, nodewidth = nodewidth, nodecol = nodecol, bordercol = bordercol)

dev.off()

####### now with prior and model selection:

hyper <- set.default.parameters(Sgenes = colnames(nem.data))

hyper$Pm <- diag(6)

print(colnames(nem.data)) # "Erk"  "Ikk2" "Jnk"  "Pi3k" "Tak1" "p38"

hyper$Pm[4:5, c(1:3,6)] <- 1

hyper$lambda = 1000
nem.res.reg <- nem(nem.data, control=hyper)

nem.res.reg <- nemModelSelection(seq(1,100, 1), nem.data, control=hyper)

plot(nem.res.reg)

tmp <- graph2adj(nem.res.reg$graph)

tmp <- transitive.reduction(tmp)

tmp2 <- adj2dnf(tmp)

count <- 0

nodelabel <- list()
nodecol <- list()
bordercol <- list()
nodeheight <- list()
nodewidth <- list()

for (i in gsub(".*=", "=", tmp2)) {
  if (length(grep(i, tmp2)) > 1) {
    count <- count + 1
    tmp2 <- gsub(i, paste("=AND", count, sep = ""), tmp2)
    tmp2 <- c(tmp2, paste("AND", count, i, sep = ""))
    nodelabel <- c(nodelabel, "")
    nodecol <- c(nodecol, "transparent")
    bordercol <- c(bordercol, "transparent")
    nodeheight <- c(nodeheight, 0)
    nodewidth <- c(nodewidth, 0)
  }
}

names(nodelabel) <- paste("AND", 1:count, sep = "")
names(nodecol) <- paste("AND", 1:count, sep = "")
names(bordercol) <- paste("AND", 1:count, sep = "")
names(nodeheight) <- paste("AND", 1:count, sep = "")
names(nodewidth) <- paste("AND", 1:count, sep = "")

pdf("publications/BNEM/BCR_nem_res_reg.pdf", height = 5, width = 3)

plotDnf(tmp2, width = 1, nodelabel = nodelabel, nodeheight = nodeheight, nodewidth = nodewidth, nodecol = nodecol, bordercol = bordercol)

dev.off()

################ other stuff:

# meld ~/comphome/publications/BNEM/manuscript/document.tex ~/comphome/Kube011_BCR_CD40_Inhibitoren/publication/pdf/manuscript/document.tex & # to check differences!
