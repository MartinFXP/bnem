## install stuff:

## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")

## BiocManager::install("CellNOptR")

## load stuff:

source("main.r")
source("low.r")
library(CellNOptR)
library(matrixStats)

## library(bnem)
## source("~/Documents/B-NEM/R/main.r"); source("~/Documents/B-NEM/R/low.r")

maxrun <- as.numeric(commandArgs(TRUE)[1])

frac <- as.numeric(commandArgs(TRUE)[2])

part <- as.numeric(commandArgs(TRUE)[3])

sd <- as.numeric(commandArgs(TRUE)[4])

n <- as.numeric(commandArgs(TRUE)[5])

## maxrun <- 2; frac <- 2; part <- 1; sd <- 0; n <- 5;

runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)

m <- 10
s <- 2
maxSize <- 2
maxStim <- 2
maxInhibit <- 1
method <- "llr"

verbose <- FALSE
draw <- FALSE

result <- array(0, c(maxrun, 5, 4), list(paste0("run", seq_len(maxrun)), c("greedy", "greedy_ia", "genetic_quick", "genetic_long", "random"), c("time", "accracy truth table", "accuracy differential effects", "score")))

for (run in runs) {
    cat(run)
    bString <- numeric(100000)
    while(length(bString) > 1000) {
        sim <- simBoolGtn(n=n, s=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, m=m, sd=sd, verbose = verbose, r = 1)
        bString <- sim$bString
        ## print(length(bString))
    }
    ## Rprof("temp.txt", line.profiling=TRUE)
    ## sim <- simBoolGtn(n=n, e=e, s=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, m=m, sd=sd, verbose = verbose, r = 1, maxcount = 1)
    ## Rprof(NULL)
    ## summaryRprof("temp.txt", lines = "show")$sampling.time
    ## head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

    ETT <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=sim$bString))

    ## plotDnf(sim$model$reacID[as.logical(reduceGraph(sim$bString, sim$model, sim$CNOlist))])

    ## plotDnf(sim$model$reacID[as.logical(reduceGraph(rand, sim$model, sim$CNOlist))])

    ## source("~/Documents/B-NEM/R/low.r"); source("~/Documents/B-NEM/R/main.r")
    start <- as.numeric(Sys.time())
    res0 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw, absorpII = FALSE)
    result[run, 1, 1] <- as.numeric(Sys.time()) - start

    ETT0 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res0$bString))

    result[run, 1, 2] <- sum(ETT0 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT0)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    tp <- (sum(sim$ERS == 1 & ERS == 1)+sum(sim$ERS == -1 & ERS == -1))
    fn <- sum(abs(sim$ERS) == 1 & sim$ERS != ERS)
    tn <- sum(sim$ERS == 0 & ERS == 0)
    fp <- sum(sim$ERS == 0 & sim$ERS != ERS)

    result[run, 1, 3] <- (tp/(tp+fn)+tn/(tn+fp))/2 # sum(ERS == sim$ERS)/length(ERS)

    result[run, 1, 4] <- min(res0$scores[[1]])

    start <- as.numeric(Sys.time())
    res1 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 2, 1] <- as.numeric(Sys.time()) - start

    ETT1 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res1$bString))

    result[run, 2, 2] <- sum(ETT1 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT1)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    tp <- (sum(sim$ERS == 1 & ERS == 1)+sum(sim$ERS == -1 & ERS == -1))
    fn <- sum(abs(sim$ERS) == 1 & sim$ERS != ERS)
    tn <- sum(sim$ERS == 0 & ERS == 0)
    fp <- sum(sim$ERS == 0 & sim$ERS != ERS)

    result[run, 2, 3] <- (tp/(tp+fn)+tn/(tn+fp))/2 # sum(ERS == sim$ERS)/length(ERS)

    result[run, 2, 4] <- min(res1$scores[[1]])

    maxTime <- result[run, 1, 1]

    start <- as.numeric(Sys.time())
    res2 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 3, 1] <- as.numeric(Sys.time()) - start

    ETT2 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res2$bString))

    result[run, 3, 2] <- sum(ETT2 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT2)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    tp <- (sum(sim$ERS == 1 & ERS == 1)+sum(sim$ERS == -1 & ERS == -1))
    fn <- sum(abs(sim$ERS) == 1 & sim$ERS != ERS)
    tn <- sum(sim$ERS == 0 & ERS == 0)
    fp <- sum(sim$ERS == 0 & sim$ERS != ERS)

    result[run, 3, 3] <- (tp/(tp+fn)+tn/(tn+fp))/2 # sum(ERS == sim$ERS)/length(ERS)

    result[run, 3, 4] <- min(res2$scores)

    maxTime <- result[run, 1, 1]*10

    start <- as.numeric(Sys.time())
    res3 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
    result[run, 4, 1] <- as.numeric(Sys.time()) - start

    ETT3 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res3$bString))

    result[run, 4, 2] <- sum(ETT3 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT3)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    result[run, 4, 3] <- sum(ERS == sim$ERS)/length(ERS)

    result[run, 4, 4] <- min(res3$scores)

    start <- as.numeric(Sys.time())
    rand <- sample(c(0,1), length(sim$model$reacID), replace = TRUE)
    result[run, 5, 1] <- as.numeric(Sys.time()) - start

    ETT4 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=rand))

    result[run, 5, 2] <- sum(ETT4 == ETT)/length(ETT)

    ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT4)
    ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]

    tp <- (sum(sim$ERS == 1 & ERS == 1)+sum(sim$ERS == -1 & ERS == -1))
    fn <- sum(abs(sim$ERS) == 1 & sim$ERS != ERS)
    tn <- sum(sim$ERS == 0 & ERS == 0)
    fp <- sum(sim$ERS == 0 & sim$ERS != ERS)

    result[run, 5, 3] <- (tp/(tp+fn)+tn/(tn+fp))/2 # sum(ERS == sim$ERS)/length(ERS)

    result[run, 5, 4] <- scoreDnf(rand, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method)

    ## result[run,,]

}

save(result, file = paste("bnem/bnem_sim", maxrun, frac, part, n, s, sd, ".rda", sep = "_"))

stop()

## leo:

module load r/3.5.1

module load curl/7.53.1

module load gmp/6.1.2

## euler:

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

##general:

ram=1000

rm error.txt

rm output.txt

rm .RData

## bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R --silent --no-save --args '2' < bnem_sim.r"

frac=20

for i in {1..20}; do
    #if [ ! -f /cluster/work/bewi/members/mpirkl/mnem_sim_results/${i}_${j}.rda ]; then
	bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R --silent --no-save --args '100' '${frac}' '${i}' '0' '10' < bnem_app.r"
    #fi
done

## plot sim:

path <- "~/Mount/Leo/"

n <- 10
s <- 2
sd <- 0
maxrun <- 100
frac <- 20
for (part in seq_len(frac)) {
    runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
    file <- paste("bnem/bnem_sim", maxrun, frac, part, n, s, sd, ".rda", sep = "_")
    if (!file.exists(paste0(path, file))) { next() }
    if (part == 1) {
        load(paste0(path, file))
        result2 <- result
    } else {
        load(paste0(path, file))
        result2[runs,,] <- result[runs,,]
    }
    cat(part)
}
result <- result2

pdf(paste("bnem_sim", n, e, s, sd, ".pdf", sep = "_"), width = 10, height = 6)
par(mfrow=c(1,4))
boxplot(result[,1:5,1], col = 2:4, ylab = "seconds", main = "Running time")
boxplot(result[,1:5,3], col = 2:4, ylab = "fraction of correct predictions", main = "Accuracy of expected differential effects")
boxplot(result[,1:5,2], col = 2:4, ylab = "fraction of correct predictions", main = "Accuracy of truth tables")
boxplot(result[,1:5,4], col = 2:4, ylab = "score between 0 and 1", main = "Score", ylim = c(0.7,1))
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

greedy0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "greedy")

greedy1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", initBstring = rep(1, length(model$reacID)))

ga0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "genetic", elitism = 1, inversion = 50)

ga1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "llr", search = "genetic", initBstring = rep(1, length(model$reacID)))

which.min(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)))

print(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)), 22)

plotDnf(ga1$graph)

plot.bnem <- function(x, ...) {
    plotDnf(res$graph)
}

bnemBs <- function(fc, x = 10, f = 0.5, replace = TRUE, independent = TRUE, startString = NULL, ...) {
    accum <- NULL
    for (i in seq_len(x)) {
        cat(i)
        fcsub <- fc[sample(seq_len(nrow(fc)), ceiling(nrow(fc)*f), replace = replace), ]
        if (is.null(startString)) {
            tmp <- bnem(fc = fcsub, ...)
        } else {
            if (!is(startString, "matrix")) {
                startString <- t(as.matrix(startString))
            }
            score <- 1
            for (j in seq_len(nrow(startString))) {
                tmp0 <- bnem(initBstring = startString[j, ], fc = fcsub, ...)
                if (min(unlist(tmp0$scores)) < score) {
                    tmp <- tmp0
                }
            }
        }

        accum <- c(accum, tmp$graph)
    }

    graph <- names(table(accum))
    freq <- table(accum)/x

    accum <- list(graph = graph, freq = freq)

    class(accum) <- "bnembs"

    return(accum)
}

initBstring = rbind(rep(0, length(model$reacID)),
                    rep(1, length(model$reacID)))

## initBstring <- initBstring[1, ]

bsres <- bnemBs(fc = fc, 1000, f = 1, CNOlist = CNOlist, model = model, method = "llr", search = "greedy", startString = initBstring, verbose = 0)

plot.bnembs <- function(x, scale = 3, shift = 0.1, cut = 0.5, dec = 2, ...) {
    graph <- x$graph
    freq <- x$freq
    graph <- graph[which(freq >= cut)]
    freq <- freq[which(freq >= cut)]
    freq2 <- NULL
    for (i in seq_len(length(graph))) {
        tmp <- rep(freq[i], length(unlist(strsplit(graph[i], "\\+"))))
        if (length(tmp) > 1) {
            freq2 <- c(freq2, tmp[1])
        }
        freq2 <- c(freq2, tmp)
    }
    freq2 <- as.vector(freq2)
    freq2 <- round(freq2, dec)
    plotDnf(graph, edgewidth = freq*scale+shift, edgelabel = freq2, ...)
}

plot(bsres, cut = 0.5, dec = 2)





