X11.options(type="Xlib")

##### old

## load("Kube011_BCR_CD40_Inhibitoren/vonMartin/BCRCD40MycCombat.RData")

## data <- BCRCD40MycCombat[, -1]

## rownames(data) <- BCRCD40MycCombat[, 1]

##### new

data <- data.combat

#### 

colnames(data) <- gsub("DMSO", "Ctrl", colnames(data))
colnames(data) <- gsub("KO", "Neg", colnames(data))
colnames(data) <- gsub("SB203580", "p38", colnames(data))
colnames(data) <- gsub("SP600125", "Jnk", colnames(data))
colnames(data) <- gsub("Ly294002", "Pi3k", colnames(data))
colnames(data) <- gsub("U0126", "Erk1xErk2", colnames(data))
colnames(data) <- gsub("IKK2", "Ikk2", colnames(data))
colnames(data) <- gsub("10058F4", "Myc", colnames(data))
colnames(data) <- gsub("TAK1", "Tak1", colnames(data))
colnames(data) <- gsub("\\+", "x", colnames(data))

data <- data[, order(colnames(data))]

batches <- paste("Batch", 1:3, sep = "")

runs <- paste("Exp", 1:2, sep = "")

stimuli <- c("BCR", "CD40")

colnames(data) <- gsub("Erk1xErk2", "Erk", colnames(data))

## colnames(data) <-gsub("Batch1", "Batch1_Exp1", gsub("Batch2", "Batch2_Exp1", gsub("Batch3", "Batch3_Exp1", colnames(data))))

## colnames(data) <-gsub("Batch4", "Batch1_Exp2", gsub("Batch5", "Batch2_Exp2", gsub("Batch6", "Batch3_Exp2", colnames(data))))

inhibitors <- c("Erk", "Ikk2", "Jnk", "p38", "Myc", "Pi3k", "Tak1")

####################### limma:

NEMlist <- list()

NEMlist$batches <- batches
NEMlist$runs <- runs
NEMlist$sep = "_"

NEMlist$exprs <- data

library(limma)
source("Boutros10.svn/method/scripts/cnopt.mod.R")
library(CellNOptR)

design <- makeDesignFull(data, stimuli, c("Vivit", "JSH", inhibitors), batches, runs, method = "raw")

for (i in 1:ncol(design)) {
  print(colnames(design)[i])
  print(colnames(data)[which(design[, i]==1)])
}

inhibitors <- c("Vivit", "JSH", inhibitors)

###################### test for ebayes prior:

## fit <- lmFit(data, design)

## pcnt <- seq(0.001,0.1, 0.001)

## de.genes <- numeric(length(pcnt))

## for (j in 1:length(pcnt)) {
##   i <- colnames(design)[2]
##   contDiff <- paste(i, "-Ctrl", sep = "")
##   contrast.matrix <- makeContrasts(contDiff, levels=design)
##   fit2 <- contrasts.fit(fit, contrast.matrix)
##   fit2 <- eBayes(fit2, proportion = pcnt[j]/100)
##   Bs <- topTable(fit2, n = nrow(data))$B
##   post.odds <- exp(Bs)/(exp(Bs)+1)
##   de.genes[j] <- sum(post.odds > 0.95)
## }

## plot(de.genes)

## lines(lowess(de.genes), col = "red")

## abline(h=median(de.genes, na.rm = T), col = "green")

########### summarize probe sets:

library(annotate)
library(hgu133plus2.db)

## data.backup <- data

data <- data.backup

rownames(data) <- mget(rownames(data), hgu133plus2SYMBOL)

data2 <- matrix(0, length(unique(rownames(data))), ncol(data))

count <- 0

for (i in unique(rownames(data))) {
  if (sum(rownames(data) %in% i) > 1) {
    count <- count + 1
    data2[count, ] <- apply(data[which(rownames(data) %in% i), ], 2, median)
  } else {
    data2[count, ] <- data[which(rownames(data) %in% i), ]
  }
}

colnames(data2) <- colnames(data)
rownames(data2) <- unique(rownames(data))
data <- data2
     
dim(data)

###################### actual foldchanges:

fit <- lmFit(data, design)

coefs <- fit$coefficients

coefs[is.na(coefs)] <- 0

NEMlist$exprs <- coefs%*%t(design)

colnames(NEMlist$exprs) <- colnames(data)

NEMlist$fc <- numeric()
NEMlist$pvals <- numeric()
NEMlist$B <- numeric()

set.prior <- 0.01 # standard

stimgrep <- intersect(grep(paste(stimuli, collapse = "|"), colnames(design)), grep(paste(inhibitors, collapse = "|"), colnames(design), invert=TRUE))
inhibgrep <- intersect(grep(paste(inhibitors, collapse = "|"), colnames(design)), grep(paste(stimuli, collapse = "|"), colnames(design), invert=TRUE))

tempNames <- character()
for (i in colnames(design)[c(stimgrep, inhibgrep)]) {
  if (is.na(fit$coefficients[1, i])) { next() }
  contDiff <- paste(i, "-Ctrl", sep = "")
  contrast.matrix <- makeContrasts(contDiff, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, proportion=set.prior)
  tempNames <- c(tempNames, paste("Ctrl_vs_", i, sep = ""))
  NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(data))$logFC[order(rownames(topTable(fit2, n = nrow(data))))])
  NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(data))$adj.P.Val[order(rownames(topTable(fit2, n = nrow(data))))])
  NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(data))$B[order(rownames(topTable(fit2, n = nrow(data))))])
}
rownames(NEMlist$fc) <- rownames(topTable(fit2, n = nrow(data)))[order(rownames(topTable(fit2, n = nrow(data))))]
rownames(NEMlist$pvals) <- rownames(topTable(fit2, n = nrow(data)))[order(rownames(topTable(fit2, n = nrow(data))))]
rownames(NEMlist$B) <- rownames(topTable(fit2, n = nrow(data)))[order(rownames(topTable(fit2, n = nrow(data))))]

for (i in colnames(design)[stimgrep]) {
  for (j in colnames(design)[inhibgrep]) {
    coefA <- paste(c(unlist(strsplit(j, "_")),unlist(strsplit(i, "_"))), collapse = "_")
    if (!(coefA %in% colnames(design)) ) { next() }
    if (is.na(fit$coefficients[1, coefA]) | is.na(fit$coefficients[1, i])) { next() }
    contDiff <- paste(coefA, "-", i, sep = "")
    ##contDiff <- paste(coefA, "-", i, "-", j, "+Ctrl", sep = "") # actual silencing effect
    contrast.matrix <- makeContrasts(contDiff, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, proportion=set.prior)
    tempNames <- c(tempNames, paste(i, "_vs_", i, "_", j, sep = ""))
    NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(data))$logFC[order(rownames(topTable(fit2, n = nrow(data))))])
    NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(data))$adj.P.Val[order(rownames(topTable(fit2, n = nrow(data))))])
    NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(data))$B[order(rownames(topTable(fit2, n = nrow(data))))])
  }
}

for (i in colnames(design)[inhibgrep]) {
  for (j in colnames(design)[stimgrep]) {
    coefA <- paste(c(unlist(strsplit(i, "_")),unlist(strsplit(j, "_"))), collapse = "_")
    if (!(coefA %in% colnames(design)) ) { next() }
    if (is.na(fit$coefficients[1, coefA]) | is.na(fit$coefficients[1, i])) { next() }
    contDiff <- paste(coefA, "-", i, sep = "")
    contrast.matrix <- makeContrasts(contDiff, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, proportion=set.prior)
    tempNames <- c(tempNames, paste(i, "_vs_", i, "_", j, sep = ""))
    NEMlist$fc <- cbind(NEMlist$fc, topTable(fit2, n = nrow(data))$logFC[order(rownames(topTable(fit2, n = nrow(data))))])
    NEMlist$pvals <- cbind(NEMlist$pvals, topTable(fit2, n = nrow(data))$adj.P.Val[order(rownames(topTable(fit2, n = nrow(data))))])
    NEMlist$B <- cbind(NEMlist$B, topTable(fit2, n = nrow(data))$B[order(rownames(topTable(fit2, n = nrow(data))))])
  }
}

colnames(NEMlist$fc) <- tempNames
colnames(NEMlist$pvals) <- tempNames
colnames(NEMlist$B) <- tempNames

NEMlist$fc <- NEMlist$fc[, sort(colnames(NEMlist$fc))]
NEMlist$pvals <- NEMlist$pvals[, sort(colnames(NEMlist$pvals))]
NEMlist$B <- NEMlist$B[, sort(colnames(NEMlist$B))]

save(NEMlist, file = "Kube011_BCR_CD40_Inhibitoren/publication/NEMlist.combat.summarized.RData")

