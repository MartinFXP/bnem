if (is.na(commandArgs(TRUE)[2])) {
    run <- as.numeric(commandArgs(TRUE)[1])
    path <- "/cluster/work/bewi/members/mpirkl/"
    library(bnem)
    library(CellNOptR)
    library(matrixStats)
    load("bcr.rda")
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
    uniquesample <- rnorm(1)
    write.table(sifMatrix, file = paste0("temp", uniquesample, ".sif"), sep = "\t", row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    PKN <- readSIF(paste0("temp", uniquesample, ".sif"))
    unlink(paste0("temp", uniquesample, ".sif"))
    CNOlist <- dummyCNOlist("BCR", c("Erk", "Ikk2", "Jnk", "p38", "Pi3k", "Tak1"), 1, 3)
    model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100, verbose = TRUE)
    initBstring = rbind(rep(0, length(model$reacID)),
                        rep(1, length(model$reacID)))
    start <- Sys.time()
    bsres <- bnemBs(fc = fc, 10, f = 1, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy", startString = initBstring, verbose = 0, draw = 0)
    end <- Sys.time()
    print(end-start)
    # save(bsres, file = paste0(path, "bnem/bcr_boot_", run, ".rda"))
    stop("bcr done")
} else {
    ## load stuff:
    library(bnem)
    library(CellNOptR)
    library(matrixStats)
    ## parameters:
    maxrun <- as.numeric(commandArgs(TRUE)[1])
    frac <- as.numeric(commandArgs(TRUE)[2])
    part <- as.numeric(commandArgs(TRUE)[3])
    maxEdges <- as.numeric(commandArgs(TRUE)[4])
    maxSize <- as.numeric(commandArgs(TRUE)[5])
    s <- as.numeric(commandArgs(TRUE)[6])
    sd <- as.numeric(commandArgs(TRUE)[7])
    n <- as.numeric(commandArgs(TRUE)[8])
    m <- as.numeric(commandArgs(TRUE)[9])
    ## maxrun <- 1; frac <- 1; part <- 1; maxEdges <- 100; maxSize <- 4; s <- 8; sd <- 1; n <- 40; m <- 10;
    runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
    maxStim <- 2
    maxInhibit <- 1
    method <- "cosine"
    verbose <- FALSE
    draw <- FALSE
    methnames <- c("greedy", "greedy_ia", "greedy_cor", "genetic_quick", "genetic_long", "genetic_stall","random")
    storenames <- c("time", "accracy truth table", "accuracy differential effects", "score","tp","fp","tn","fn")
    result <- array(0, c(maxrun, length(methnames), length(storenames)), list(paste0("run", seq_len(maxrun)), methnames, storenames))
    path <- "/cluster/work/bewi/members/mpirkl/"
    for (run in runs) {
        cat(run)
        bString <- numeric(100000)
        while(length(bString) > 1000) {
            sim <- simBoolGtn(Sgenes=n, maxEdges = maxEdges, stimGenes=s, maxSize = maxSize, maxStim=maxStim, maxInhibit=maxInhibit, Egenes=m, verbose = verbose, reps = 1, frac = 0.1, layer = 3, negation = 0.33, and = 0.25)
            bString <- sim$bString
            cat(".")
        }
        ## mnem::plotDnf(sim$PKN$reacID)
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
        pos <- which(sim$fc == 1)
        neg <- which(sim$fc == -1)
        zero <- which(sim$fc == 0)
        sim$fc[pos] <- rnorm(length(pos),1,sd)
        sim$fc[neg] <- rnorm(length(neg),-1,sd)
        sim$fc[zero] <- rnorm(length(zero),0,sd)
        start <- as.numeric(Sys.time())
        res0 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw, absorpII = FALSE)
        result[run, 1, 1] <- as.numeric(Sys.time()) - start
        ETT0 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res0$bString))
        result[run, 1, 2] <- sum(ETT0 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT0)
        ERS0 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 1, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 1, 4] <- min(res0$scores[[1]])
        result[run, 1, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 1, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 1, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 1, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        start <- as.numeric(Sys.time())
        res1 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
        result[run, 2, 1] <- as.numeric(Sys.time()) - start
        ETT1 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res1$bString))
        result[run, 2, 2] <- sum(ETT1 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT1)
        ERS1 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 2, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 2, 4] <- min(res1$scores[[1]])
        result[run, 2, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 2, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 2, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 2, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        maxTime <- result[run, 2, 1]
        start <- as.numeric(Sys.time())
        res1 <- bnem(search = "greedy", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = "s", verbose = verbose, draw = draw)
        result[run, 3, 1] <- as.numeric(Sys.time()) - start
        ETT1 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res1$bString))
        result[run, 3, 2] <- sum(ETT1 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT1)
        ERS1 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 3, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 3, 4] <- min(res1$scores[[1]])
        result[run, 3, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 3, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 3, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 3, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        start <- as.numeric(Sys.time())
        res2 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw,stallGenMax=Inf)
        result[run, 4, 1] <- as.numeric(Sys.time()) - start
        ETT2 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res2$bString))
        result[run, 4, 2] <- sum(ETT2 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT2)
        ERS2 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 4, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 4, 4] <- min(res2$scores)
        result[run, 4, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 4, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 4, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 4, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        maxTime <- result[run, 2, 1]*10
        start <- as.numeric(Sys.time())
        res3 <- bnem(search = "genetic", maxTime = maxTime, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw,stallGenMax=Inf)
        result[run, 5, 1] <- as.numeric(Sys.time()) - start
        ETT3 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res3$bString))
        result[run, 5, 2] <- sum(ETT3 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT3)
        ERS3 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 5, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 5, 4] <- min(res3$scores)
        result[run, 5, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 5, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 5, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 5, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        start <- as.numeric(Sys.time())
        res4 <- bnem(search = "genetic", fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method, verbose = verbose, draw = draw)
        result[run, 6, 1] <- as.numeric(Sys.time()) - start
        ETT3 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=res4$bString))
        result[run, 6, 2] <- sum(ETT3 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT3)
        ERS3 <- ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 6, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 6, 4] <- min(res4$scores)
        result[run, 6, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 6, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 6, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 6, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
        start <- as.numeric(Sys.time())
        rand <- sample(c(0,1), length(sim$model$reacID), replace = TRUE)
        result[run, 7, 1] <- as.numeric(Sys.time()) - start
        ETT4 <- t(simulateStatesRecursive(CNOlist=sim$CNOlist, model=sim$model, bString=rand))
        result[run, 7, 2] <- sum(ETT4 == ETT)/length(ETT)
        ERS <- computeFc(CNOlist=sim$CNOlist, y = ETT4)
        ERS <- ERS[, which(colnames(ERS) %in% colnames(sim$ERS))]
        result[run, 7, 3] <- sum(ERS == sim$ERS)/length(ERS)
        result[run, 7, 4] <- scoreDnf(rand, fc = sim$fc, CNOlist = sim$CNOlist, model = sim$model, method = method)
        result[run, 7, 5] <- sum(ERS == 1 & sim$ERS == 1)+sum(ERS==-1 & sim$ERS == -1)
        result[run, 7, 6] <- sum(abs(ERS) == 1 & sim$ERS == 0)
        result[run, 7, 7] <- sum(ERS == 0 & sim$ERS == 0)
        result[run, 7, 8] <- sum(ERS == 0 & abs(sim$ERS) == 1)
    }
    save(result, file = paste0(path, paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_")))
    stop("simulation done")
}

## HPC commands to run simulations:

## scp bnem_app.r mpirkl@euler:

ram=1000
rm error.txt
rm output.txt
rm .RData
maxedges=100
maxin=4
R=/cluster/work/bewi/members/mpirkl/miniconda3/bin/R

Sgenes=25
queue=4
frac=10
Egenes=10
Stimuli=2
Noise=1

if [ ${Sgenes} == '20' ]
then
Stimuli=4
frac=50
ram=4000
fi
if [ ${Sgenes} == '25' ]
then
Stimuli=5
frac=50
ram=6000
fi
if [ ${Sgenes} == '30' ]
then
Stimuli=6
frac=100
ram=8000
queue=24
fi
if [ ${Sgenes} == '40' ]
then
Stimuli=8
frac=100
ram=16000
queue=24
fi

ram=$(expr $ram \* 2)

## maxrun frac part maxedges maxgatesize stims noise sgenes egenes

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${R} --vanilla --silent --no-save --args '100' '${frac}' '1' '${maxedges}' '${maxin}' '${Stimuli}' '${Noise}' '${Sgenes}' '${Egenes}' < bnem_app.r"

for i in $(eval echo {2..$frac}) ## {2..100}; do
do
    #if [ ! -f /cluster/work/bewi/members/mpirkl/mnem_sim_results/${i}_${j}.rda ]; then
    if [[ $i == 3 || $i == 26 || $i == 27 || $i == 29 || $i == 33 || $i == 36 || $i == 39 || $i == 40 || $i == 43 || $i == 46 || $i == 50 ]]; then
	bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${R} --vanilla --silent --no-save --args '100' '${frac}' '${i}' '${maxedges}' '${maxin}' '${Stimuli}' '${Noise}' '${Sgenes}' '${Egenes}' < bnem_app.r"
	fi
done

## bootstrap:

ram=1000
rm error.txt
rm output.txt
rm .RData

for i in {1..100}; do
        bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --vanilla --no-save --args '${i}' < bnem_app.r"
done

## plot results:

library(mnem)
path <- "Z:/" # "~/Mount/Eulershare/"
n <- 10
s <- 2
m <- 10
sd <- 1
frac <- 50
results <- list()
count <- 0
maxrun <- 100
Sgenes <- c(10,20,25,30,40)
for (n in Sgenes) {
    result2 <- NULL
    if (n==10) {s<-2;frac<-10} else if (n==20) {s<-4;frac<-50} else if (n==25) {s<-5;frac<-50} else if (n==30) {s<-6;frac<-100} else if (n==40) {s<-8;frac<-100}
    for (part in seq_len(frac)) {
        runs <- (maxrun/frac*part - maxrun/frac + 1):(maxrun/frac*part)
        file <- paste0(path,paste("bnem/bnem_sim", n, m, s, sd, maxrun, frac, part, ".rda", sep = "_"))
        if (!file.exists(file)) { cat(paste0(part,"."));next() }
        load(file)
        if (dim(result)[3] < 8) {
            print(file)
            next()
        }
        if (is.null(result2)) {
            result2 <- result
        } else {
            result2[runs,,] <- result[runs,,]
        }
    }
    result <- result2
    if (is.null(result)) { result <- results[[count]] }
    count <- count + 1
    results[[count]] <- result
}

n.meth <- dim(results[[1]])[2]
if (n.meth==7) {
    methods <- c("greedy","greedy\n(inverse absorption)","greedy\n(spearman)","genetic\n(time-limit)","genetic","genetic (stall)","random")
    cols <- c("darkred","red","#FF8888","darkgreen","green","lightgreen","grey")
} else {
    methods <- c("greedy","greedy\n(inverse absorption)","genetic\n(time-limit)","genetic","genetic (stall)","random")
    cols <- c("darkred","red","darkgreen","green","lightgreen","grey")
}

wilcox <- array(NA,c(length(results),length(methods),length(methods)),list(Sgenes=Sgenes,methods=methods,methods2=methods))
idx1 <- 5
idx2 <- 8
pw <- TRUE
for (i in 1:length(Sgenes)) {
    for (j in 1:n.meth) {
        for (k in 1:n.meth) {
            if (j==k) next()
        wilcox[i,j,k] <- wilcox.test(results[[i]][,j,idx1]/(results[[i]][,j,idx1]+results[[i]][,j,idx2]),results[[i]][,k,idx1]/(results[[i]][,k,idx1]+results[[i]][,k,idx2]),alternative="less",paired=pw)$p.value
        }
    }
}

myboxplot <- mnem:::myboxplot
box <- 1
time <- 1
idxs <- c(1,3,5)
if (box) {
    restime <- cbind(results[[idxs[1]]][,1:n.meth,1],results[[idxs[2]]][,1:n.meth,1],results[[idxs[3]]][,1:n.meth,1])
    if (time) {
        pdf("temp.pdf", width = 11, height = 3)
        laymat <- matrix(c(rep(1,50),rep(2,50),rep(3,50),rep(4,32)),1,byrow=TRUE)
    } else {
        pdf("temp.pdf", width = 12, height = 3)
        laymat <- matrix(c(rep(1,50),rep(2,50),rep(3,20)),1,byrow=TRUE)
        print(apply(restime,2,median))
        print(apply(restime,2,mean))
    }
    layout(laymat)
    v.idx <- c(n.meth+0.5,n.meth*2+0.5)
    axis.idx <- c(n.meth/2+0.5,n.meth*1.5+0.5,n.meth*2.5+0.5)
    if (time) {
        myboxplot(restime+10^-10, col = cols,border=cols,medcol="black",ylab = "seconds (log10-scale)", main = "Running time", box = box,dens=0,xaxt="n",bordercol=cols,log="y",
                  ylim=c(max(min(restime[,-(n.meth*(1:3))]),10^-1),max(restime[,-(n.meth*(1:3))])))
        abline(v=v.idx)
        axis(1,axis.idx,Sgenes[idxs])
    }
    resacc <- cbind(results[[idxs[1]]][,1:n.meth,3],results[[idxs[2]]][,1:n.meth,3],results[[idxs[3]]][,1:n.meth,3])
    resacc <- cbind(results[[idxs[1]]][,1:n.meth,5],results[[idxs[2]]][,1:n.meth,5],results[[idxs[3]]][,1:n.meth,5])/(cbind(results[[idxs[1]]][,1:n.meth,5],results[[idxs[2]]][,1:n.meth,5],results[[idxs[3]]][,1:n.meth,5])+cbind(results[[idxs[1]]][,1:n.meth,6],results[[idxs[2]]][,1:n.meth,6],results[[idxs[3]]][,1:n.meth,6]))
    myboxplot(resacc, col = cols,border=cols,medcol="black",ylab = "precision", main = "Precision of expected differential effects", box = box,dens=0,xaxt="n",bordercol=cols,ylim=c(0,1))
    abline(v=v.idx)
    axis(1,axis.idx,Sgenes[idxs])
    resll <- -cbind(results[[idxs[1]]][,1:n.meth,4],results[[idxs[2]]][,1:n.meth,4],results[[idxs[3]]][,1:n.meth,4])
    resll <- cbind(results[[idxs[1]]][,1:n.meth,5],results[[idxs[2]]][,1:n.meth,5],results[[idxs[3]]][,1:n.meth,5])/(cbind(results[[idxs[1]]][,1:n.meth,5],results[[idxs[2]]][,1:n.meth,5],results[[idxs[3]]][,1:n.meth,5])+cbind(results[[idxs[1]]][,1:n.meth,8],results[[idxs[2]]][,1:n.meth,8],results[[idxs[3]]][,1:n.meth,8]))
    myboxplot(resll, col = cols,border=cols,medcol="black",ylab = "recall", main = "Recall of expected differential effects", box = box,dens=0,xaxt="n",bordercol=cols,ylim=c(0,1))#, ylim = c(0,1),dens=0,bordercol=cols)
    abline(v=v.idx)
    axis(1,axis.idx,Sgenes[idxs])
    op <- par(mar = rep(0, 4))
    plot(1:100,col="transparent",yaxt="n",xaxt="n",bty="n",xlab="",ylab="")
    legend('top',legend=methods,col=cols,fill=cols,border=FALSE,box.lwd=0,box.col="transparent",y.intersp=2)
    par(op)
    dev.off()
}

## analyze BCR:

data(bcr) ## load data
bcr <- processDataBCR() ## alternative
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
devtools::load_all("../../")
greedy0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy")
greedy1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy", initBstring = rep(1, length(model$reacID)))
ga0 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "genetic",draw=0)
ga1 <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "genetic", initBstring = rep(1, length(model$reacID)),draw=0)
which.min(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)))
print(c(min(greedy0$scores[[1]]), min(greedy1$scores[[1]]), min(ga0$scores), min(ga1$scores)), 22)
plotDnf(greedy1$graph)
initBstring = rbind(rep(0, length(model$reacID)),
                    rep(1, length(model$reacID)))

if (min(greedy0$scores[[1]]) < min(greedy0$scores[[1]])) {
    bcrbest <- greedy0
} else {
    bcrbest <- greedy1
}

## reduce further:
cbind(bcrbest$bString,model$reacID)
bestString <- bcrbest$bString
bestString[12] <- 0
bestString[9] <- 1
par(mfrow=c(1,2))
plotDnf(model$reacID[as.logical(bcrbest$bString)])
plotDnf(model$reacID[as.logical(bestString)])
scoreDnf(bcrbest$bString,CNOlist,fc,model=model)
scoreDnf(bestString,CNOlist,fc,model=model)

## read and plot bootstrap results:
bsfull <- NULL
for (i in 1:100) {
    file <- paste0("~/Mount/Eulershare/bnem_old/bcr_boot_", i, ".rda")
    if (!file.exists(file)) { cat(i); next() }
    load(file)
    bsfull <- c(bsfull, bsres$x)
    cat(".")
}
## check the result:
initString <- numeric(length(model$reacID))
initString[model$reacID %in% names(table(bsfull))[table(bsfull)>0.5*1000]] <- 1
greedyBS <- bnem(fc = fc, CNOlist = CNOlist, model = model, method = "cosine", search = "greedy",initBstring = initString)
bsfull <- list(x = toupper(bsfull), n = 1000)
class(bsfull) <- "bnemBs"
dnf2016 <- c("BCR=TAK1","BCR=PI3K","PI3K=JNK","TAK1=ERK","TAK1=IKK2","PI3K=IKK2","JNK=p38","PI3K+IKK2=p38")
library(bnem)
pdf("temp.pdf", width = 12, height = 6)
par(mfrow=c(1,3))
plotDnf(dnf2016, nodeshape = list(BCR = "diamond"),edgelwd=3)
plotDnf(model$reacID[as.logical(bestString)], nodeshape = list(BCR = "diamond"),edgelwd=3)
plot(bsfull, cut = 0.1, dec = 2, ci = 0, nodeshape = list(BCR = "diamond"),
     fontsize = 14,labelcol='blue')
dev.off()

## compare old (2016) and new bootstrap:
dnf2016 <- c("BCR=Pi3k","BCR=Tak1","Pi3k=Ikk2","Pi3k=Jnk","Tak1=Ikk2","Tak1=Erk","Jnk=p38","Pi3k+Ikk2=p38")
dnfBoot <- c("BCR=Pi3k","BCR=Tak1","Pi3k=Erk","Pi3k=Jnk","Tak1=Ikk2","Tak1=Erk","Tak1+Ikk2+Jnk=p38")
string2016 <- numeric(length(model$reacID))
string2016[model$reacID %in% dnf2016] <- 1
stringBoot <- numeric(length(model$reacID))
stringBoot[model$reacID %in% dnfBoot] <- 1
scoreDnf(string2016,CNOlist,fc,model=model)
scoreDnf(stringBoot,CNOlist,fc,model=model)
scoreDnf(string2016,CNOlist,fc,model=model,method="s")
scoreDnf(stringBoot,CNOlist,fc,model=model,method="s")
par(mfrow=c(2,2))
plotDnf(model$reacID[as.logical(string2016)])
plotDnf(model$reacID[as.logical(stringBoot)])












combiNames2 <- paste(rownames(CNOlist@stimuli), combiInhibit, sep = "_")[which(rowSums(cbind(CNOlist@stimuli, CNOlist@inhibitors)[grepStimsKds, ]) != 0)]

combiNames <- NULL
for (i in grepStimsKds) {
    combiNames <- c(combiNames, paste(names(which(stimsKdsCbind[i, ] >= 1)), collapse = "_"))
}

combiNames2 <- combiNames0[grepStimsKds]

cbind(combiNames, combiNames2)

all(combiNames == combiNames2)
