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

## bsres <- bnemBs(fc = fc, 10, f = 1, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", startString = initBstring, verbose = 0)

## save(bsres, file = paste0("bnem/bcr_boot_", run, ".rda"))

## stop()

## load stuff:

print(version)

source("main.r")
source("low.r")
library(CellNOptR)
library(matrixStats)

## library(bnem)
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

## maxrun <- 10; frac <- 1; part <- 1; maxEdges <- 1000; maxSize <- 2; s <- 2; sd <- 1; n <- 10; m <- 10

runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)

maxStim <- 2
maxInhibit <- 1
method <- "llr"

verbose <- FALSE
draw <- FALSE

result <- array(0, c(maxrun, 5, 4), list(paste0("run", seq_len(maxrun)), c("greedy", "greedy_ia", "genetic_quick", "genetic_long", "random"), c("time", "accracy truth table", "accuracy differential effects", "score")))

for (run in runs) {

    ## run <- 1

    cat(run)
    bString <- numeric(100000)
    while(length(bString) > 1000) {
        cat(".")
        sim <- simBoolGtn(Sgenes=n, maxEdges = maxEdges, stimGenes=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, Egenes=m, sd=sd, verbose = verbose, reps = 1, frac = 0.1)
        bString <- sim$bString
    }
    ## plotDnf(sim$PKN$reacID)
    ## Rprof("temp.txt", line.profiling=TRUE)
    ## sim <- simBoolGtn(n=n, e=e, s=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, m=m, sd=sd, verbose = verbose, r = 1, maxcount = 1)
    ## Rprof(NULL)
    ## summaryRprof("temp.txt", lines = "show")$sampling.time
    ## head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

    ETT <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=sim$bString))

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

frac=100

## maxrun frac part maxedges maxgatesize stims noise sgenes egenes

for i in {1..1}; do
    #if [ ! -f /cluster/work/bewi/members/mpirkl/mnem_sim_results/${i}_${j}.rda ]; then
	bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '100' '${frac}' '${i}' '1000' '2' '2' '1' '10' '10' < bnem_app.r"
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

path <- "~/Mount/Leo/"

n <- 10
s <- 5
sd <- 1
m <- 2
maxrun <- 100
frac <- 100
for (part in seq_len(frac)) {
    runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
    file <- paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_")
    if (!file.exists(paste0(path, file))) { cat(part);next() }
    if (part == 1) {
        load(paste0(path, file))
        result2 <- result
    } else {
        load(paste0(path, file))
        result2[runs,,] <- result[runs,,]
    }
}
result <- result2

pdf("temp.pdf", width = 10, height = 10)
par(mfrow=c(1,4))
boxplot(result[,1:5,1], col = 2:4, ylab = "seconds", main = "Running time")
boxplot(result[,1:5,3], col = 2:4, ylab = "fraction of correct predictions", main = "Accuracy of expected differential effects")
boxplot(result[,1:5,2], col = 2:4, ylab = "fraction of correct predictions", main = "Accuracy of truth tables")
boxplot(result[,1:5,4], col = 2:4, ylab = "score between 0 and 1", main = "Score")#, ylim = c(0,1))
dev.off()

## paper fig:

pdf("temp.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
boxplot(result[,2:4,1], col = 2:5, ylab = "seconds", main = "running time", xaxt = "n")
axis(1, 1:4, c("Greedy", "Gen_s", "Gen_l", "rand"))
boxplot(1-result[,2:4,4], col = 2:5, ylab = "probability", main = "normalized likelihood", xaxt = "n")#, ylim = c(0,0.3))
axis(1, 1:4, c("Greedy", "Gen_s", "Gen_l", "rand"))
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

greedy0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "greedy")

greedy1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", initBstring = rep(1, length(model$reacID)))

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
