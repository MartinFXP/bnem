# main functions:
source("Boutros10.svn/method/scripts/computeFc.R")
source("Boutros10.svn/method/scripts/computeFcII.R")
source("Boutros10.svn/method/scripts/computeFcII2.R")
source("Boutros10.svn/method/scripts/computeFcIII.R")
source("Boutros10.svn/method/scripts/computeFcDisc.R")
source("Boutros10.svn/method/scripts/computeFcAsp.R")
source("Boutros10.svn/method/scripts/computeSm.R")
source("Boutros10.svn/method/scripts/getNemFit.R")
source("Boutros10.svn/method/scripts/gaBinaryNemT1.R")
source("Boutros10.svn/method/scripts/computeScoreNemT1.R")
source("Boutros10.svn/method/scripts/localSearch.R")
source("Boutros10.svn/method/scripts/preClingo.R") # is probably not being used at the moment
source("Boutros10.svn/method/scripts/gibbsSamp.R")
source("Boutros10.svn/method/scripts/gaRecon.R")
source("Boutros10.svn/method/scripts/simulateStates.R")
source("Boutros10.svn/method/scripts/simulateStatesRecursive.R")
source("Boutros10.svn/method/scripts/simulateStatesDAG.R")
source("Boutros10.svn/method/scripts/rSimulatorT1.R")
source("Boutros10.svn/method/scripts/reduceGraph.R")
source("Boutros10.svn/method/scripts/expandNEM.R")
source("Boutros10.svn/method/scripts/checkNEMlist.R")
source("Boutros10.svn/method/scripts/checkCNOlist.R")
source("Boutros10.svn/method/scripts/simulateStatesRecursiveIII.R")
source("Boutros10.svn/method/scripts/simulateStatesRecursiveII.R")
source("Boutros10.svn/method/scripts/hbNEM.R")
source("Boutros10.svn/method/scripts/partitionGraph.R")
source("Boutros10.svn/method/scripts/exSearch.R")
source("Boutros10.svn/method/scripts/simulateStatesRecursiveAdd.R")
source("Boutros10.svn/method/scripts/crossTalk.R")
source("Boutros10.svn/method/scripts/sbnem.R")
source("Boutros10.svn/method/scripts/snem.R")
source("Boutros10.svn/method/scripts/triples.R")
source("Boutros10.svn/method/scripts/isDag.R")
source("Boutros10.svn/method/scripts/dnf2adj.R")
source("Boutros10.svn/method/scripts/removeCycles.R")
source("Boutros10.svn/method/scripts/absorption.R")
source("Boutros10.svn/method/scripts/BP.R")
source("Boutros10.svn/method/scripts/mc.R")
source("Boutros10.svn/method/scripts/resolveTrue.R")
source("Boutros10.svn/method/scripts/getHierarchy.R")
source("Boutros10.svn/method/scripts/checkMethod.R")
source("Boutros10.svn/method/scripts/drawLocal.R")
source("Boutros10.svn/method/scripts/exportVars.R")
source("Boutros10.svn/method/scripts/plotDnf.R")
source("Boutros10.svn/method/scripts/simulateDnf.R")
source("Boutros10.svn/method/scripts/resBNEM.R")
source("Boutros10.svn/method/scripts/addEdge.R")
source("Boutros10.svn/method/scripts/deleteEdge.R")
source("Boutros10.svn/method/scripts/addSignal.R")
source("Boutros10.svn/method/scripts/deleteSignal.R")
source("Boutros10.svn/method/scripts/dataRed.R")

#preprocessing:
source("Boutros10.svn/method/scripts/optScore.R")
source("Boutros10.svn/method/scripts/cnoFromData.R")
source("Boutros10.svn/method/scripts/egeneSel.R")
source("Boutros10.svn/method/scripts/preproData.R")
source("Boutros10.svn/method/scripts/addDumpnodes.R")
source("Boutros10.svn/method/scripts/dummyCNOlist.R")
source("Boutros10.svn/method/scripts/norm.R")
source("Boutros10.svn/method/scripts/makeDesign.R")
source("Boutros10.svn/method/scripts/makeDesign2.R")
source("Boutros10.svn/method/scripts/makeDesignFull.R")
source("Boutros10.svn/method/scripts/makeSif.R")
##source("Boutros10.svn/src/classes/KnockdownSet.R")

#result analysis:
source("Boutros10.svn/method/scripts/heatmapRaw.R")
source("Boutros10.svn/method/scripts/plotBinary.R")
source("Boutros10.svn/method/scripts/validateGraph.R")
source("Boutros10.svn/method/scripts/residualErrorNem.R")
source("Boutros10.svn/method/scripts/drawScores.R")
source("Boutros10.svn/method/scripts/myGsea.R")
source("Boutros10.svn/method/scripts/heatmapOP.R")
source("Boutros10.svn/method/scripts/convertGraph.R")
source("Boutros10.svn/method/scripts/transClose.R")
source("Boutros10.svn/method/scripts/transRed.R")

# get internal original functions: (v 1.4)

simulatorT0 <- get("simulatorT0", en = asNamespace("CellNOptR"))
addPriorKnowledge <- get("addPriorKnowledge", en = asNamespace("CellNOptR"))
cutModel <- get("cutModel", en = asNamespace("CellNOptR"))

createCube <- function(n=3, m=n) {
  if (m > n) { m <- n }
  n2 <- FALSE
  m2 <- FALSE
  if (!(round(n/2) == n/2)) {
    n2 <- TRUE
    n <- n + 1
  }
  if (!(round(m/2) == m/2)) {
    m2 <- TRUE
    m <- m + 1
  }
  e <- matrix(0, m, n)
  e[1, ] <- sample(c(0,1), n, replace = T)
  for (i in (1:(m/2))*2) {
    e[i, ] <- 1 - e[i-1, ]
    if (i >= m) { break() }
    for (j in 1:i) {
      e[i+1, ((j-1)*(n/i)+1):(j*(n/i))] <- e[j, ((j-1)*(n/i)+1):(j*(n/i))]
    }
  }
  if (n2) {
    e <- e[, -1]
  }
  if (m2) {
    e <- e[-nrow(e), ]
  }
  return(e)
}

                                  
gsub2 <- function(pattern, replace, x) {
  p <- list()
  y <- x
  for (i in 1:length(pattern)) {
    p[[i]] <- regexp(pattern[i], x)
  }
  
}
.int2dec <- function(x){
    return(sum(x*2^(rev(seq(x))-1)))
}

my.sim <- function(x,y, m = "c", c1 = 0.25, c2 = 0.75, s1 = 1, s2 = 1, s3 = 1) { 
  if (m %in% "c") {
    return((x%*%y)/(((x%*%x)^0.5)*((y%*%y)^0.5)))
  }
  if (m %in% "p") {
    return(((x-mean(x))%*%(y-mean(y)))/((((x-mean(x))%*%(x-mean(x)))^0.5)*(((y-mean(y))%*%(y-mean(y)))^0.5)))
  }
  if (m %in% "s") {
    return(((rank(x)-mean(rank(x)))%*%(rank(y)-mean(rank(y))))/((((rank(x)-mean(rank(x)))%*%(rank(x)-mean(rank(x))))^0.5)*(((rank(y)-mean(rank(y)))%*%(rank(y)-mean(rank(y))))^0.5)))
  }
  if (m %in% "k") {
    cp <- 0
    dp <- 0
    x.r <- rank(x)
    y.r <- rank(y)
    n <- length(x)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if ((x.r[i] < x.r[j] & y.r[i] < y.r[j]) | (x.r[i] > x.r[j] & y.r[i] > y.r[j])) {
          cp <- cp + 1
        }
        if ((x.r[i] < x.r[j] & y.r[i] > y.r[j]) | (x.r[i] > x.r[j] & y.r[i] < y.r[j])) {
          dp <- dp + 1
        }
      }
    }
    tau <- (cp - dp)/(0.5*n*(n-1))
    return(tau)
  }
  if (m %in% "mad") {
    S0 <- as.matrix(1 - abs(x))
    Spos <- x
    Spos[which(Spos == -1)] <- 0
    Sneg <- x
    Sneg[which(Sneg == 1)] <- 0
    E0 <- y
    E0[which(y == c1)] <- s1
    E0[which(y > c1 & y < s2)] <- s2*s1
    E0[which(y >= s2)] <- -s3*s1
    Epos <- y
    Epos[which(y == c1)] <- -s3*s1*s2
    Epos[which(y > c1 & y < s2)] <- s2
    Epos[which(y >= s2)] <- 1
    Epos[which(y < c1 & y > -s2)] <- -s3*s2
    Epos[which(y <= -s2)] <- -s3
    Eneg <- y
    Eneg[which(y == c1)] <- -s3*s1*s2
    Eneg[which(y < c1 & y > -s2)] <- s2
    Eneg[which(y <= -s2)] <- 1
    Eneg[which(y > c1 & y < -s2)] <- -s3*s2
    Eneg[which(y >= -s2)] <- -s3
    EposI <- -y
    EposI[which(y == c1)] <- -s3*s1*s2
    EposI[which(y > c1 & y < s2)] <- s2
    EposI[which(y >= s2)] <- 1
    EposI[which(y < c1 & y > -s2)] <- -s3*s2
    EposI[which(y <= -s2)] <- -s3
    EnegI <- -y
    EnegI[which(y == c1)] <- -s3*s1*s2
    EnegI[which(y < c1 & y > -s2)] <- s2
    EnegI[which(y <= -s2)] <- 1
    EnegI[which(y > c1 & y < -s2)] <- -s3*s2
    EnegI[which(y >= -s2)] <- -s3
    ES0 <- E0%*%S0
    ESpos <- Epos%*%Spos
    ESneg <- Eneg%*%Sneg
    ESposI <- EposI%*%Spos
    ESnegI <- EnegI%*%Sneg
    return(max(ES0 + ESpos + ESneg, ES0 + ESposI + ESnegI)/length(x))
  }
}

X.prod <- function(x,y) {

  if (length(x) != 3 | length(y) != 3) {
     return("vectors have to be of length 3")
   }

  z <- numeric(3)
  z[1] <- x[2]*y[3] - x[3]*y[2]
  z[2] <- x[3]*y[1] - x[1]*y[3]
  z[3] <- x[1]*y[2] - x[2]*y[1]

  return(z)
}

fisher.cor.test <- function(x,y, method = "p") {
  r <- cor(x,y, method = method)
  z <- atanh(r)*(length(x) - 3)^0.5
  return(2*pnorm(-abs(z)))
}

perm.test <- function(x,y, size = 100, alternative = "two.sided", Z = NULL, method = "p") {
  samples <- numeric(size)
  n <- length(x)
  if (is.null(Z)) {
    for (i in 1:size) {
      x2 <- x
      y2 <- y[sample(1:n, n)]
      samples[i] <- cor(x2, y2, method = method)
    }
    if (alternative %in% "two.sided") {
      p.value <- sum(abs(samples) >= abs(cor(x,y, method = method)))/size
    }
    if (alternative %in% "greater") {
      p.value <- sum(samples >= abs(cor(x,y, method = method)))/size
    }
    if (alternative %in% "less") {
      p.value <- sum(abs(samples) <= abs(cor(x,y, method = method)))/size
    }
  } else {
    for (i in 1:size) {
      x2 <- x
      y2 <- y[sample(1:n, n)]
        z2 <- Z[sample(1:n, n)]
        samples[i] <- pcor.test(x,y2,z2, method = method)$estimate
    }
    if (alternative %in% "two.sided") {
      p.value <- sum(abs(samples) >= abs(pcor.test(x,y,Z, method = method)$estimate))/size
    }
    if (alternative %in% "greater") {
      p.value <- sum(samples >= abs(pcor.test(x,y,Z, method = method)$estimate))/size
    }
    if (alternative %in% "less") {
      p.value <- sum(abs(samples) <= abs(pcor.test(x,y,Z, method = method)$estimate))/size
    }
  }
  return(p.value)
}

par.cor <- function(x,y,Z, method = "p") {
 r12 <- cor(x,y, method = method)
 r13 <- cor(x,Z, method = method)
 r23 <- cor(y,Z, method = method)
 r <- (r12 - r13*r23)/(((1-r13)^0.5)*((1-r23)^0.5))
 return(r)
}


whiskplot <- function(x, width = 1, length = 1, lines = FALSE, ylim = NULL, range = sd, values = median, gap = 0.1, ...) {
  x.values <- apply(x, 2, values)
  x.mean <- apply(x, 2, mean)
  x.range <- apply(x, 2, range)
  if (lines == FALSE) {
    if (is.null(ylim)) {
      ylim <- c(min(x.mean)-max(x.range), max(x.mean)+max(x.range))
    }
    plot(x.values, ylim = ylim, ...)
  } else {
    lines(x.values, ...)
  }
  for (i in 1:length(x.mean)) {
    
    segments(y0=x.values[i] + gap*2*max(x.range), y1=(x.mean[i] + x.range[i]*width), x0=i, x1=i, ...)
    
    segments(y0=(x.mean[i] - x.range[i]*width), y1=x.values[i] - gap*2*max(x.range), x0=i, x1=i, ...)
    
    ##segments(y0=(x.mean[i] - x.range[i]*width), y1=(x.mean[i] + x.range[i]*width), x0=i, x1=i, ...)
    
    segments(x0=(i - length), x1=(i + length), y0=(x.mean[i] + x.range[i]*width), y1=(x.mean[i] + x.range[i]*width), ...)
    segments(x0=(i - length), x1=(i + length), y0=(x.mean[i] - x.range[i]*width), y1=(x.mean[i] - x.range[i]*width), ...)
   
  }
}

## ToothGrowth$dose.cat <- factor(ToothGrowth$dose, labels=paste("d", 1:3, sep=""))
## df <- with(ToothGrowth , aggregate(len, list(supp=supp, dose=dose.cat), mean))
## df$se <- with(ToothGrowth , aggregate(len, list(supp=supp, dose=dose.cat), 
##               function(x) sd(x)/sqrt(10)))[,3]

## opar <- theme_update(panel.grid.major = theme_blank(),
##                      panel.grid.minor = theme_blank(),
##                      panel.background = theme_rect(colour = "black"))

## xgap <- position_dodge(0.2)
## gp <- ggplot(df, aes(x=dose, y=x, colour=supp, group=supp))
## gp + geom_line(aes(linetype=supp), size=.6, position=xgap) + 
##      geom_point(aes(shape=supp), size=3, position=xgap) + 
##      geom_errorbar(aes(ymax=x+se, ymin=x-se), width=.1, position=xgap)
## theme_set(opar)

                   
