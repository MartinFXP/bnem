# main functions:
source("github/trunk/method/computeFc.R")
source("github/trunk/method/computeFcII.R")
source("github/trunk/method/computeFcII2.R")
source("github/trunk/method/computeFcIII.R")
source("github/trunk/method/computeFcDisc.R")
source("github/trunk/method/computeFcAsp.R")
source("github/trunk/method/computeSm.R")
source("github/trunk/method/getNemFit.R")
source("github/trunk/method/gaBinaryNemT1.R")
source("github/trunk/method/computeScoreNemT1.R")
source("github/trunk/method/localSearch.R")
source("github/trunk/method/preClingo.R") # is probably not being used at the moment
source("github/trunk/method/gibbsSamp.R")
source("github/trunk/method/gaRecon.R")
source("github/trunk/method/simulateStates.R")
source("github/trunk/method/simulateStatesRecursive.R")
source("github/trunk/method/simulateStatesDAG.R")
source("github/trunk/method/rSimulatorT1.R")
source("github/trunk/method/reduceGraph.R")
source("github/trunk/method/expandNEM.R")
source("github/trunk/method/checkNEMlist.R")
source("github/trunk/method/checkCNOlist.R")
source("github/trunk/method/simulateStatesRecursiveIII.R")
source("github/trunk/method/simulateStatesRecursiveII.R")
source("github/trunk/method/hbNEM.R")
source("github/trunk/method/partitionGraph.R")
source("github/trunk/method/exSearch.R")
source("github/trunk/method/simulateStatesRecursiveAdd.R")
source("github/trunk/method/crossTalk.R")
source("github/trunk/method/sbnem.R")
source("github/trunk/method/snem.R")
source("github/trunk/method/triples.R")
source("github/trunk/method/isDag.R")
source("github/trunk/method/dnf2adj.R")
source("github/trunk/method/adj2dnf.R")
source("github/trunk/method/removeCycles.R")
source("github/trunk/method/absorption.R")
source("github/trunk/method/BP.R")
source("github/trunk/method/mc.R")
source("github/trunk/method/resolveTrue.R")
source("github/trunk/method/getHierarchy.R")
source("github/trunk/method/checkMethod.R")
source("github/trunk/method/drawLocal.R")
source("github/trunk/method/exportVars.R")
source("github/trunk/method/plotDnf.R")
source("github/trunk/method/compressDnf.R")
source("github/trunk/method/simulateDnf.R")
source("github/trunk/method/resBNEM.R")
source("github/trunk/method/addEdge.R")
source("github/trunk/method/deleteEdge.R")
source("github/trunk/method/addSignal.R")
source("github/trunk/method/deleteSignal.R")
source("github/trunk/method/dataRed.R")
source("github/trunk/method/compressDnf.R")

#preprocessing:
source("github/trunk/method/optScore.R")
source("github/trunk/method/cnoFromData.R")
source("github/trunk/method/egeneSel.R")
source("github/trunk/method/preproData.R")
source("github/trunk/method/addDumpnodes.R")
source("github/trunk/method/dummyCNOlist.R")
source("github/trunk/method/norm.R")
source("github/trunk/method/makeDesign.R")
source("github/trunk/method/makeDesign2.R")
source("github/trunk/method/makeDesignFull.R")
source("github/trunk/method/makeSif.R")
##source("Boutros10.svn/src/classes/KnockdownSet.R")

#result analysis:
source("github/trunk/method/heatmapRaw.R")
source("github/trunk/method/plotBinary.R")
source("github/trunk/method/validateGraph.R")
source("github/trunk/method/residualErrorNem.R")
source("github/trunk/method/drawScores.R")
source("github/trunk/method/myGsea.R")
source("github/trunk/method/heatmapOP.R")
source("github/trunk/method/convertGraph.R")
source("github/trunk/method/transClose.R")
source("github/trunk/method/transRed.R")
source("github/trunk/method/smoothMatrix.R")
source("github/trunk/method/randomDnf.R")

# get internal original functions: (v 1.4)

simulatorT0 <- get("simulatorT0", en = asNamespace("CellNOptR"))
addPriorKnowledge <- get("addPriorKnowledge", en = asNamespace("CellNOptR"))
cutModel <- get("cutModel", en = asNamespace("CellNOptR"))

gene2protein <- function(genes, strict = FALSE) {
  if (strict) {
    gene2prot <- cbind(
                   c("^APC$", "^ATF2$", "^BIRC2$", "^BIRC3$", "^CASP4$", "^CASP8$", "^CFLAR$", "^CHUK$", "^CTNNB1$", "^DKK1$", "^DKK4$", "^FLASH$", "^IKBKB$", "^IKBKG$", "^JUN$", "^MAP2K1$", "^MAP3K14$", "^MAP3K7$", "^MAPK8$", "^PIK3CA$", "^RBCK1$", "^RELA$", "^RIPK1$", "^RIPK3$", "^RNF31$", "^SHARPIN$", "^TAB2$", "^TCF4$", "^TCF7L2$", "^TNFRSF10A$", "^TNFRSF10B$", "^TNFRSF1A$", "^TNIK$", "^TRAF2$", "^USP2$", "^WLS$", "^WNT11$", "^WNT5A$", "^TNFa$", "^TRAIL"),
                   c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                   )
  } else {
    gene2prot <- cbind(
                   c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                   c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                   )
  }
  for (i in 1:nrow(gene2prot)) {
    genes <- gsub(gene2prot[i, 1], gene2prot[i, 2], genes)
  }
  return(genes)
}

protein2gene <- function(proteins, strict = FALSE) {
  if (strict) {
    gene2prot <- cbind(
                   c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                   c("^Apc$", "^Atf2$", "^cIap1$", "^cIap2$", "^Casp4$", "^Casp8$", "^c-Flip$", "^Ikka$", "^Beta-Cat.$", "^Dkk1$", "^Dkk4$", "^Casp8ap2$", "^Ikkb$", "^Nemo$", "^cJun$", "^Mekk$", "^Nik$", "^Tak1$", "^Jnk$", "^Pi3k$", "^Hoil1$", "^RelA$", "^Rip1$", "^Rip3$", "^Hoip$", "^Sharpin$", "^Tab2$", "^fake$", "^Tcf4$", "^Dr4$", "^Dr5$", "^Tnfr1$", "^Tnik$", "^Traf2$", "^Usp2$", "^Evi$", "^Wnt11$", "^Wnt5A$", "^Tnfa$", "^Trail$")
                   )
  } else {
    gene2prot <- cbind(
                   c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8", "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH", "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                   c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa", "Trail")
                   )
  }
  for (i in 1:nrow(gene2prot)) {
    proteins <- gsub(gene2prot[i, 2], gene2prot[i, 1], proteins)
  }
  return(proteins)
}
    
myCN2bioCN <- function(x, stimuli, inhibitors) {
  y <- gsub("_vs_", ") vs (", x)
  for (i in inhibitors) {
    y <- gsub(i, paste(i, "\\-", sep = ""), y)
  }
  for (i in stimuli) {
    y <- gsub(i, paste(i, "\\+", sep = ""), y)
  }
  y <- gsub("Ctrl", "control", paste("(", gsub("_", ",", y), ")", sep = ""))
  return(y)
}

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

bin2dec <- function(x) {
  exp2 <- 2^c((length(x)-1):0)
  y <- exp2%*%x
  return(y)
}

dec2bin <- function(x) {
  if (x == 0) {
    y <- 0
  } else {
    tmp <- rev(as.integer(intToBits(x)))
    y <- tmp[which(tmp != 0)[1]:length(tmp)] # this is faster than my dec2bin
  }
  return(y)
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

                   


    
