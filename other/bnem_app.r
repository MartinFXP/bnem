## install stuff:

## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")

## BiocManager::install("CellNOptR")

## bootstraps:

## run <- as.numeric(commandArgs(TRUE)[1])

## source("main.r")
## source("low.r")
## library(CellNOptR)
## library(matrixStats)

## load("bcr.rda")

## fc <- bcr$fc

## sifMatrix <- rbind(c("BCR", 1, "Pi3k"),
##                    c("BCR", 1, "Tak1"),
##                    c("Tak1", 1, "Erk"),
##                    c("Pi3k", 1, "Erk"),
##                    c("Tak1", 1, "p38"),
##                    c("Pi3k", 1, "p38"),
##                    c("Ikk2", 1, "p38"),
##                    c("Jnk", 1, "p38"),
##                    c("Pi3k", 1, "Ikk2"),
##                    c("Tak1", 1, "Ikk2"),
##                    c("p38", 1, "Ikk2"),
##                    c("Jnk", 1, "Ikk2"),
##                    c("Pi3k", 1, "Jnk"),
##                    c("Tak1", 1, "Jnk"),
##                    c("Ikk2", 1, "Jnk"),
##                    c("p38", 1, "Jnk"))
## uniquesample <- rnorm(1)
## write.table(sifMatrix, file = paste0("temp", uniquesample, ".sif"), sep = "\t", row.names = FALSE,
##             col.names = FALSE,
##             quote = FALSE)
## PKN <- readSIF(paste0("temp", uniquesample, ".sif"))
## unlink(paste0("temp", uniquesample, ".sif"))

## CNOlist <- dummyCNOlist("BCR", c("Erk", "Ikk2", "Jnk", "p38", "Pi3k", "Tak1"), 1, 3)

## model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100, verbose = TRUE)

## initBstring = rbind(rep(0, length(model$reacID)),
##                     rep(1, length(model$reacID)))

## ## initBstring <- initBstring[1, ]

## nemfc <- bcr$fc[, c(1,2,6,8,9,10)]
## colnames(nemfc) <- gsub(".*_", "", colnames(nemfc))
## nemfc[which(abs(nemfc) >= log2(1.5))] <- 1
## nemfc[which(nemfc != 1)] <- 0
## nemres <- mnem:::mynem(nemfc, method = "disc")

## bsres <- bnemBs(fc = fc, 10, f = 1, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy", startString = initBstring, verbose = 0)

## save(bsres, file = paste0("bnem/bcr_boot_", run, ".rda"))

## stop()

## load stuff:

print(version)

source("main.r")
source("low.r")
library(CellNOptR)
library(matrixStats)

## source("~/Documents/B-NEM/R/main.r"); source("~/Documents/B-NEM/R/low.r")

maxrun <- as.numeric(commandArgs(TRUE)[1])

frac <- as.numeric(commandArgs(TRUE)[2])

part <- as.numeric(commandArgs(TRUE)[3])

maxEdges <- as.numeric(commandArgs(TRUE)[4])

maxSize <- as.numeric(commandArgs(TRUE)[5])

s <- as.numeric(commandArgs(TRUE)[6])

sd <- as.numeric(commandArgs(TRUE)[7])

n <- as.numeric(commandArgs(TRUE)[8])

m <- as.numeric(commandArgs(TRUE)[9])

## maxrun <- 10; frac <- 1; part <- 1; maxEdges <- 100; maxSize <- 2; s <- 6; sd <- 0.5; n <- 30; m <- 10

runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)

maxStim <- 2
maxInhibit <- 1
method <- "cosine"

verbose <- FALSE
draw <- FALSE

result <- array(0, c(maxrun, 5, 4), list(paste0("run", seq_len(maxrun)), c("greedy", "greedy_ia", "genetic_quick", "genetic_long", "random"), c("time", "accracy truth table", "accuracy differential effects", "score")))

for (run in runs) {

    ## run <- 1

    ## source("~/Documents/B-NEM/R/low.r"); source("~/Documents/B-NEM/R/main.r")
    cat(run)
    bString <- numeric(100000)
    while(length(bString) > 1000) {
        sim <- simBoolGtn(Sgenes=n, maxEdges = maxEdges, stimGenes=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, Egenes=m, sd=sd, verbose = verbose, reps = 3, frac = 0.1, layer = 3, negation = 0, and = 0.25)
        bString <- sim$bString
        cat(".")
    }

    double <- NULL
    for (i in colnames(sim$CNOlist@stimuli)) {
        for (j in colnames(sim$CNOlist@stimuli)) {
            if (i == j) { next() }
            double <- c(double, paste(sort(c(i,j)), collapse = "_"))
        }
    }
    double <- unique(double)

    sim$fc <- sim$fc[, grep(paste(c(paste0("^", colnames(sim$CNOlist@stimuli), "_vs"),
                                    paste0("^", double, "_vs"),
                                    paste0("Ctrl_vs_", colnames(sim$CNOlist@stimuli), "$"),
                                    paste0("Ctrl_vs_", double, "$")
                                    ),
                                  collapse = "|"), colnames(sim$fc))]
    sim$ERS <- sim$ERS[, unique(colnames(sim$fc))]

    ETT <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=sim$bString))

    ## par(mfrow=c(1,2)); plotDnf(sim$model$reacID[as.logical(sim$bString)]); plotDnf(sim$PKN$reacID)

    ## ## runtime:
    ## source("~/Documents/B-NEM/R/low.r"); source("~/Documents/B-NEM/R/main.r")
    ## Rprof("temp.txt", line.profiling=TRUE)
    ## res0 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw, absorpII = FALSE, maxSteps = 1)
    ## Rprof(NULL)
    ## summaryRprof("temp.txt", lines = "show")$sampling.time
    ## head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

    ## source("~/Documents/B-NEM/R/low.r"); source("~/Documents/B-NEM/R/main.r")
    start <- as.numeric(Sys.time())
    res0 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw, absorpII = FALSE)
    result[run, 1, 1] <- as.numeric(Sys.time()) - start

    ETT0 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res0$bString))

    result[run, 1, 2] <- sum(ETT0 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT0)
    ERS0 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 1, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 1, 4] <- min(res0$scores[[1]])

    start <- as.numeric(Sys.time())
    res1 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 2, 1] <- as.numeric(Sys.time()) - start

    ETT1 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res1$bString))

    result[run, 2, 2] <- sum(ETT1 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT1)
    ERS1 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 2, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 2, 4] <- min(res1$scores[[1]])

    maxTime <- result[run, 2, 1]

    start <- as.numeric(Sys.time())
    res2 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 3, 1] <- as.numeric(Sys.time()) - start

    ETT2 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res2$bString))

    result[run, 3, 2] <- sum(ETT2 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT2)
    ERS2 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 3, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 3, 4] <- min(res2$scores)

    maxTime <- result[run, 2, 1]*10

    start <- as.numeric(Sys.time())
    res3 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 4, 1] <- as.numeric(Sys.time()) - start

    ETT3 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res3$bString))

    result[run, 4, 2] <- sum(ETT3 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT3)
    ERS3 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 4, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 4, 4] <- min(res3$scores)

    start <- as.numeric(Sys.time())
    rand <- sample(c(0,1), length(sim$model$reacID), replace = TRUE)
    result[run, 5, 1] <- as.numeric(Sys.time()) - start

    ETT4 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=rand))

    result[run, 5, 2] <- sum(ETT4 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT4)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 5, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 5, 4] <- scoreDnf(rand, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method)

    ## result[1,,]; par(mfrow=c(1,4)); plotDnf(sim$model$reacID[as.logical(res1$bString)]); plotDnf(sim$model$reacID[as.logical(sim$bString)]); plotDnf(sim$model$reacID[as.logical(res2$bString)]); plotDnf(sim$model$reacID[as.logical(res3$bString)]);

    ## result[run,,]

}

save(result, file = paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_"))

stop()

## leo:

module load r/3.5.1

module load curl/7.53.1

module load gmp/6.1.2

## euler: (use R/bin/R instead of R)

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

##general:

ram=1000

rm error.txt

rm output.txt

rm .RData

## bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R --silent --no-save --args '2' < bnem_app.r"

queue=4

frac=100

Sgenes=30
Egenes=10
Stimuli=6
Noise=1

## maxrun frac part maxedges maxgatesize stims noise sgenes egenes

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '100' '${frac}' '1' '100' '2' '${Stimuli}' '${Noise}' '${Sgenes}' '${Egenes}' < bnem_app.r"

for i in $( eval echo {2..$frac} ) ## {2..100}; do
do
    #if [ ! -f /cluster/work/bewi/members/mpirkl/mnem_sim_results/${i}_${j}.rda ]; then
	bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '100' '${frac}' '${i}' '100' '2' '${Stimuli}' '${Noise}' '${Sgenes}' '${Egenes}' < bnem_app.r"
    #fi
done

## bootstrap:

ram=1000

rm error.txt

rm output.txt

rm .RData

for i in {1..100}; do
        bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R --silent --no-save --args '${i}' < bnem_app.r"
done

## plot sim:

path <- "~/Mount/Leo/" # path <- "~/Mount/Euler/"

n <- 30
m <- 10
s <- 6
sd <- 1
frac <- 100

result2 <- NULL
maxrun <- 100
for (part in seq_len(frac)) {
    runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
    file <- paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_")
    if (!file.exists(paste0(path, file))) { cat(part);next() }
    if (is.null(result2)) {
        load(paste0(path, file))
        result2 <- result
    } else {
        load(paste0(path, file))
        result2[runs,,] <- result[runs,,]
    }
}
result <- result2

box <- 1
source("mnem/R/mnems_low.r")
pdf("temp.pdf", width = 10, height = 10)
par(mfrow=c(1,4))
cols <- rgb(c(1,0,1,0,1), c(0,1,0,1,0), c(1,1,0,0,0), 0.75)
myboxplot(result[,1:5,1], col = cols, ylab = "seconds", main = "Running time", box = box)
myboxplot(result[,1:5,3], col = cols, ylab = "fraction of correct predictions", main = "Accuracy of expected differential effects", box = box)
myboxplot(result[,1:5,2], col = cols, ylab = "fraction of correct predictions", main = "Accuracy of truth tables", box = box)
myboxplot(result[,1:5,4], col = cols, ylab = "log likelihood + constant", main = "likelihood", box = box)#, ylim = c(0,1))
dev.off()

## paper fig(s):

pdf("temp.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
boxplot(result[,1:5,1], col = 2:5, ylab = "seconds", main = "running time", xaxt = "n")
axis(1, 1:5, dimnames(result)[[2]])
boxplot(-result[,1:5,4], col = 2:5, ylab = "log likelihood + constant", main = "likelihood", xaxt = "n")#, ylim = c(0,0.3))
axis(1, 1:5, dimnames(result)[[2]])
dev.off()

## more in one:

path <- "~/Mount/Leo/" # path <- "~/Mount/Euler/"

ns <- 10
m <- 10
s <- 2
sds <- c(0.5, 1, 2)
sdl <- length(sds)
fracs <- c(100,10,10)

results <- list()
for (n in ns) {
    for (sd in sds) {
        result2 <- NULL
        maxrun <- 100
        frac <- fracs[which(sds == sd)]
        for (part in seq_len(frac)) {
            runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
            file <- paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_")
            if (!file.exists(paste0(path, file))) { cat(part); next() }
            if (is.null(result2)) {
                load(paste0(path, file))
                result2 <- result
            } else {
                load(paste0(path, file))
                result2[runs,,] <- result[runs,,]
            }
        }
        results[[which(sds == sd)]] <- result2
    }

    tmp0 <- tmp1 <- NULL
    for (i in 1:length(results)) {
        tmp0 <- cbind(tmp0, results[[i]][,2:4,1])
        tmp1 <- cbind(tmp1, results[[i]][,2:4,4])
    }
}

cols <- rgb(c(1,0,0),c(0,1,0),c(0,0,1), 0.75)
pdf("temp.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
boxplot(tmp0, col = cols, ylab = "seconds", main = "running time", xaxt = "n")
axis(1, 1:(sdl*3), rep(c("Greedy", "Gen_S", "Gen_L"), sdl))
abline(v=c(3.5, 6.5)[1:(length(sds)-1)], lty = 2, col = rgb(0,0,0,0.75))
axis(1, c(2,5,8)[1:length(sds)], sds, line = 2, tick = 0)
boxplot(-tmp1, col = cols, ylab = "log likelihood + constant", main = "network score", xaxt = "n")
axis(1, 1:(sdl*3), rep(c("Greedy", "Gen_S", "Gen_L"), sdl))
abline(v=c(3.5, 6.5)[1:(length(sds)-1)], lty = 2, col = rgb(0,0,0,0.75))
axis(1, c(2,5,8)[1:length(sds)], sds, line = 2, tick = 0)
dev.off()

## analyze BCR:

data(bcr)

fc <- bcr$fc

sifMatrix <- rbind(c("BCR", 1, "Pi3k"),
                   c("BCR", 1, "Tak1"),
                   c("Tak1", 1, "Erk"),
                   c("Pi3k", 1, "Erk"),
                   c("Tak1", 1, "p38"),
                   c("Pi3k", 1, "p38"),
                   c("Ikk2", 1, "p38"),
                   c("Jnk", 1, "p38"),
                   c("Pi3k", 1, "Ikk2"),
                   c("Tak1", 1, "Ikk2"),
                   c("p38", 1, "Ikk2"),
                   c("Jnk", 1, "Ikk2"),
                   c("Pi3k", 1, "Jnk"),
                   c("Tak1", 1, "Jnk"),
                   c("Ikk2", 1, "Jnk"),
                   c("p38", 1, "Jnk"))
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
PKN <- readSIF("temp.sif")
unlink('temp.sif')

CNOlist <- dummyCNOlist("BCR", c("Erk", "Ikk2", "Jnk", "p38", "Pi3k", "Tak1"), 1, 3)

model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100, verbose = TRUE)

source("~/Documents/B-NEM/R/main.r"); source("~/Documents/B-NEM/R/low.r")

greedy0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy")

greedy1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy", initBstring = rep(1, length(model$reacID)))

ga0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "genetic")

ga1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "genetic", initBstring = rep(1, length(model$reacID)))

greedyM <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", initBstring = ga0$bString)

which.min(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)))

print(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)), 22)

plotDnf(greedy0$graph)

initBstring = rbind(rep(0, length(model$reacID)),
                    rep(1, length(model$reacID)))

## initBstring <- initBstring[1, ]

## takes long use hpc

bsres <- bnemBs(fc = fc, 10, f = 0.5, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", startString = initBstring, verbose = 0)

## read hpc results:

bsfull <- NULL

for (i in 1:100) {
    if (!file.exists(paste0("~/Mount/Leo/bnem/bcr_boot_", i, ".rda"))) { cat(i); next() }
    load(paste0("~/Mount/Leo/bnem/bcr_boot_", i, ".rda"))
    bsfull <- c(bsfull, bsres)
}

bsfull <- list(x = bsfull, n = 1000)
class(bsfull) <- "bnembs"

pdf("temp.pdf", width = 8, height = 8)
plot(bsfull, cut = 0.5, dec = 2, ci = 0, nodeshape = list(BCR = "diamond"))
dev.off()

bcr2 <- processDataBCR(path = "~/Downloads/celfiles/", combsign = 0)
















combiNames2 <- paste(rownames(CNOlist@stimuli), combiInhibit, sep = "_")[which(rowSums(cbind(CNOlist@stimuli, CNOlist@inhibitors)[grepStimsKds, ]) != 0)]

combiNames <- NULL
for (i in grepStimsKds) {
    combiNames <- c(combiNames, paste(names(which(stimsKdsCbind[i, ] >= 1)), collapse = "_"))
}

combiNames2 <- combiNames0[grepStimsKds]

cbind(combiNames, combiNames2)

all(combiNames == combiNames2)
