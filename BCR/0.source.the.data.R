X11.options(type="Xlib")

source("Boutros10.svn/method/scripts/cnopt.mod.R")

## new:

load("Kube011_BCR_CD40_Inhibitoren/expression_set.rdat")

head(expr.set)

data <- exprs(expr.set)

library(sva)

colnames(data) <- gsub("Batch$", "Batch3", gsub("IKK2", "Ikk2", gsub("LY", "Ly", gsub("Ko", "KO", colnames(data)))))

batches <- gsub("Batch$", "Batch3", gsub("\\.", "", gsub(".*Batch", "Batch", colnames(data))))

stimuli <- c("BCR", "CD40")

inhibitors <- unique(sort(gsub("_.*", "", colnames(data))))

inhibitors <- inhibitors[-2]

both <- c(stimuli, inhibitors)

design <- matrix(0, length(batches), length(both))

for (i in 1:length(both)) {
  tmp <- numeric(length(batches))
  tmp[grep(both[i], colnames(data))] <- 1
  design[, i] <- tmp
}

design <- cbind(1, design)

design <- rep(1, length(batches))

stimuli <- c("BCR", "CD40")

inhibitors <- c("TAK1", "Ly294002", "Ikk2", "SP600125", "SB203580", "U0126", "10058F4", "Tak1", "JSH", "Vivit")

batches2 <- c("Batch1", "Batch2", "Batch3")

runs <- c("Exp1", "Exp2")

colnames(data) <- gsub("Batch6", "Batch3_Exp2", gsub("Batch5", "Batch2_Exp2", gsub("Batch4", "Batch1_Exp2", gsub("Batch3", "Batch3_Exp1", gsub("Batch2", "Batch2_Exp1", gsub("Batch1", "Batch1_Exp1", colnames(data)))))))

data <- data[, order(colnames(data))]

source("Boutros10.svn/method/scripts/cnopt.mod.R")

design <- makeDesignFull(data, stimuli, inhibitors, batches2, runs = runs, method = "raw")

design <- design[, -grep("Ctrl|Batch", colnames(design))]

for (i in 1:ncol(design)) {
  print(colnames(design)[i])
  print(colnames(data)[which(design[, i] == 1)])
}

genes.mean <- apply(data, 1, mean)

data <- data[order(genes.mean, decreasing = T)[1:10000], ]

genes.sd <- apply(data, 1, sd)

data <- data[order(genes.sd, decreasing = T)[1:1000], ]

data.combat <- ComBat(data, batches, rep(1, ncol(data)))#, design)

d <- numeric(ncol(data))

colseq <- order(d, decreasing = T)
#colseq <- 1:ncol(data)

par(ask=T)
for (i in colseq) {
  plot(data[, i], data.combat[, i], main = colnames(data)[i])
  abline(0, 1, col = "red")
  d[i] <- dist(rbind(data[, i], data.combat[, i]))
}
