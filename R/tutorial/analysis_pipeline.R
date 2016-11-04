###### this is a documented guide to allow anyone to do a B-NEM analysis on their own dataset:

X11.options(type="Xlib")

source(".../cnopt.mod.R") # load all functions; adjust your working directory or the paths in cnopt.mod.R # source("github/trunk/method/cnopt.mod.R")

library(CellNOptR) # CNO package required

## with your own data you obviously skip the next step and continue with "NEMlist <- list()".

data <- matrix(rnorm(100*5), 100, 5) # this should be your data
rownames(data) <- 1:100
colnames(data) <- c("Ctrl_vs_S", paste("S_vs_S_", LETTERS[1:4], sep = "")) # typical NEM setup with one stimulation S and four single knock-downs A-D

# NOTE: stimulations and inhibited S-genes can overlap. E.g. if you have one S-gene X, which is stimulated in one experiment and inhibited in another. Rename X to Xi and add the node Xs. Later you must add only the edge Xs -> Xi c("Xs", 1, "Xi") in your PKN. Xs is not directly connected to anything else! This way you do have two different S-genes Xs and Xi in your model, but in practice it works as if you only have X and it can be stimulated in one and inhibited in another experiment. Note, that the inhibition of X overrules the stimulation in this scenario.

NEMlist <- list()
NEMlist$exprs <- NULL
NEMlist$fc <- data # these are your foldchanges, e.g. limma contrasts (Smyth et al.) with featers (genes/mRNA) as rows and contrasts as of different conditions as columns.

## The data has to be anntated correctly. For the contrast "S1 - control" the corresponding sample/column is named
## "Ctrl_vs_S1" (= stimulation effect, if S1 is a stimulation). For the contrast "S1_Sk - S1" the
## sample/column has to be named "S1_vs_S1_Sk" (=correspondign silencing effect of Sk during S1 stimulation). "S1_Sk" is the condition with S1 stimulated and Sk inhibited. Look at the toy_example.R for a more detailed description

stimuli <- c("S")
inhibitors <- LETTERS[1:4]
CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1)

Sgenes <- c(stimuli, inhibitors)

## if you do want to do full reconstruction:
sifMatrix <- numeric()
for (i in Sgenes) {
  for (j in Sgenes) {
    if (i %in% j) { next() }
    sifMatrix <- rbind(sifMatrix, c(i, 1, j))
    sifMatrix <- rbind(sifMatrix, c(i, -1, j)) # if you want negative edges
  }
}
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

## if you want to use prior knowledge:
PKN <- readSIF(".../yourPKN.sif") # sif is tab-delimited. Every line has one interaction. A 1 B (=A activates B). A -1 B (=A inhibits B).

model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100)

## now comes the inference:

initBstring <- rep(0, length(model$reacID)) # start with the empty network
initBstring <- absorption(rep(1, length(model$reacID))) # or fully connected OR-gates only.

parallel <- 8 # list(c(4,16,8,2), c("machine1", "machine2", "machine3", "machine4")) # parallel calculation of edge improvements

# verbose = FALSE can be set for all search algorithms

# try greedy search:

locRun <- localSearch(
         CNOlist=CNOlist,
         NEMlist=NEMlist,
         model=model,
         parallel=parallel,
         initSeed=initBstring,
         draw = TRUE # FALSE does not draw the network evolution and can be faster
         )

resString <- locRun$bStrings[1, ]

# or genetic algorithm:

gaRun <- gaBinaryNemT1(
           parallel = parallel,
           CNOlist=CNOlist,
           NEMlist = NEMlist,
           model=model,
           initBstring=initBstring,
           popSize = 100,
           stallGenMax = 10,
           graph = TRUE # FALSE does not draw the network evolution and can be faster
           )

resString <- gaRun$bString

# exhaustive search (only recommended for very "small" search spaces; < 20 hyper-edges):

exRun <- exSearch(
           parallel = parallel,
           CNOlist=CNOlist,
           NEMlist = NEMlist,
           model=model
           )

resString <- exRun$bString

# lets look at the data and how well it fits the resolved network:
# with experience you can use this to identify missing edges in your prior network. obviously not in the toy example since the GTN is a subnetwork of the extended prior.

geneLists <- list()
par(ask=T)
for (i in 1:ncol(CNOlist@signals[[1]])) {
  geneLists[[i]] <- validateGraph(CNOlist, NEMlist, model = model, bString = resString, Sgene = i, Egenes = 1000, cexRow = 0.8, soft = T, cexCol = 0.7, xrot = 45, disc = 0, Colv = T, Rowv = T, dendrogram = "both", bordercol = "grey", aspect = "iso", sub = "")
  dev.print("temp.pdf", device = pdf, width = 40, height = 10) # take a closer look
}; names(geneLists) <- colnames(CNOlist@signals[[1]]); par(ask=F)

# have fun with your own analysis

# one important parameter, which can be trained before the final optimization starts is zeta. This controls the sparseness of the solution and is set to 10^-10 by default. You can change that with "sizeFac=a" with any number a in the optimization functions localSearch or gaBinaryNemT1.

# training zeta:

results <- list()

targets <- list()

count <- 0

zetas <- c(1, rev(seq(10^-5,1,0.2)^2), 0)

plot(c(rev(seq(10^-5,1,length.out=100)^2), 0), type = "l")

plot(zetas, type = "b")

trainruns <- 2

for (i in 1:trainruns) {

  NEMlist2 <- NEMlist

  targets[[i]] <- sample(1:nrow(NEMlist2$fc), floor(0.5*nrow(NEMlist2$fc)))

  NEMlist2$fc <- NEMlist2$fc[targets[[i]], ]

  for (j in zetas) {

    par(ask=F, mfrow=c(1,1))
    count <- count + 1
    results[[count]] <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist2, model=model, parallel=parallel, initSeed = rep(0, length(model$reacID)), sizeFac = j)

  }

}

## save(results, CNOlist, model, NEMlist, targets, file = "temp.RData")

load("temp.RData")

rep <- matrix(0, length(zetas), trainruns)
pred <- matrix(0, length(zetas), trainruns)
graph.size <- matrix(0, length(zetas), trainruns)
node.num <- matrix(0, length(zetas), trainruns)

count <- 0

for (i in 1:trainruns) {

  NEMlist2 <- NEMlist
  NEMlist2$fc <- NEMlist2$fc[targets[[i]], ]
  NEMlist3 <- NEMlist
  NEMlist3$fc <- NEMlist3$fc[-targets[[i]], ]

  for (j in 1:length(zetas)) {

    count <- count + 1
    
    pred[j, i] <- computeScoreNemT1(CNOlist, model = model, results[[count]]$bString, NEMlist = NEMlist3, sizeFac = 0)
    rep[j, i] <- computeScoreNemT1(CNOlist, model = model, results[[count]]$bString, NEMlist = NEMlist2, sizeFac = zetas[j])
    graph.tmp <- model$reacID[as.logical(results[[count]]$bString)]
    if (length(graph.tmp) > 0) {
      nodes.tmp <- length(unique(gsub("!", "", unlist(strsplit(unlist(strsplit(graph.tmp, "=")), "\\+")))))
    } else {
      nodes.tmp <- 0
    }
    graph.size[j, i] <- length(unlist(strsplit(graph.tmp, "\\+")))
    node.num[j, i] <- nodes.tmp

  }

}

if (min(rowMeans(rep)) < min(rowMeans(pred))) {
  ymin <- min(rowMeans(rep))
  rep2 <- rowMeans(rep) - min(rowMeans(rep))
  pred2 <- rowMeans(pred) - min(rowMeans(rep))
} else {
  ymin <- min(rowMeans(pred))
  rep2 <- rowMeans(rep) - min(rowMeans(pred))
  pred2 <- rowMeans(pred) - min(rowMeans(pred))
}
if (max(rep2) > max(pred2)) {
  ymax <- max(rowMeans(rep))
  pred2 <- pred2/max(rep2)
  rep2 <- rep2/max(rep2)
} else {
  ymax <- max(rowMeans(pred))
  rep2 <- rep2/max(pred)
  pred2 <- pred2/max(pred2)
}

graph.size2 <- rowMeans(graph.size)/length(model$reacID)

node.num2 <- rowMeans(node.num)/length(c(stimuli, inhibitors))

plot(1 - pred2, type = "b", xaxt = "n", yaxt = "n", xlab = "zetas", ylab = "test set scores", pch = "T", ylim = c(0,1), sub = "test set scores (T), train set scores (L), number of connected nodes (N) in %, graph size (G) in %")
axis(1, 1:7, c(round(zetas, 2)[1:5], zetas[6:7]))
axis(2, seq(1,100,length.out=10)/100, round(seq(ymin, ymax, length.out=10), 2))
lines(1 - rep2, type = "b", pch = "L")
lines(graph.size2, type = "b", pch = "G")
lines(node.num2, type = "b", pch = "N")
axis(4, seq(1,100,length.out=10)/100, seq(1, 100, length.out=10))
abline(v=which.min(pred2), col = "red", lty = 3)
abline(h=min(pred2), col = "blue", lty = 3)
abline(h=graph.size2[which.min(pred2)], col = "blue", lty = 3)
abline(h=node.num2[which.min(pred2)], col = "blue", lty = 3)
