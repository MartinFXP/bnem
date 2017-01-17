#' @noRd
computeSm <-
function (CNOlist,
                       NEMlist,
                       parameters,
                       method="standard"
                       ) {
  CompMatCont <- NEMlist$fc
  fitScore <- -1 # match of high degree
  errorScore <- parameters$scoring[3] # mismatch of high degree
  fitMult <- parameters$scoring[2] # multiplicator for match of low degree
  errorMult <- parameters$scoring[2] # multiplicator for mismatch of low degree
  zeroMult <- parameters$scoring[1]
  
  CompMat <- CompMatCont
  Epos <- CompMat
  Eneg <- CompMat
  E0 <- CompMat
  EposI <- CompMat*(-1)
  EnegI <- CompMat*(-1)

  beta <- parameters$cutOffs[2]
  alpha <- parameters$cutOffs[1]

  wrongPos <- which(Epos < -beta)
  rightPosII <- which(Epos > alpha & Epos < beta)
  rightPosI <- which(Epos >= beta)
  zerosPosI <- which(Epos <= alpha & Epos >= -alpha)
  zerosPosII <- which(Epos < -alpha & Epos > -beta)
  
  wrongNeg <- which(Eneg > beta)
  rightNegII <- zerosPosII # which(Eneg < -alpha & Eneg > -beta)
  rightNegI <- which(Eneg <= -beta)
  zerosNegI <- zerosPosI # which(Eneg >= -alpha & Eneg <= alpha)
  zerosNegII <- rightPosII # which(Eneg > alpha & Eneg < beta)
  
  right0I <- which(abs(E0) <= alpha)
  right0II <- which(abs(E0) > alpha & abs(E0) < beta)
  wrong0 <- which(abs(E0) >= beta)
  
  if ("mLL" %in% method | "cp" %in% method) {
    
    E0 <- E0*0
    E0[right0I] <- 1
    
    Epos <- Epos*0
    Epos[rightPosI] <- 1
    
    Eneg <- Eneg*0
    Eneg[rightNegI] <- 1
    
    EposI <- EposI*0
    EposI[rightPosII] <- 1
    
    EnegI <- EnegI*0
    EnegI[rightNegII] <- 1
    
  } else {

    E0[wrong0] <- errorScore*zeroMult
    E0[right0II] <- fitScore*zeroMult*fitMult
    E0[right0I] <- fitScore*zeroMult
    
    Epos[zerosPosI] <- errorScore*errorMult*zeroMult
    Epos[zerosPosII] <- errorScore*errorMult
    Epos[wrongPos] <- errorScore
    Epos[rightPosI] <- fitScore
    Epos[rightPosII] <- fitScore*fitMult
    
    Eneg[zerosNegI] <- -errorScore*errorMult*zeroMult
    Eneg[zerosNegII] <- -errorScore*errorMult
    Eneg[wrongNeg] <- -errorScore
    Eneg[rightNegI] <- -fitScore
    Eneg[rightNegII] <- -fitScore*fitMult
    
    EposI[zerosNegI] <- errorScore*errorMult*zeroMult
    EposI[zerosNegII] <- errorScore*errorMult
    EposI[wrongNeg] <- errorScore
    EposI[rightNegI] <- fitScore
    EposI[rightNegII] <- fitScore*fitMult
    
    EnegI[zerosPosI] <- -errorScore*errorMult*zeroMult
    EnegI[zerosPosII] <- -errorScore*errorMult
    EnegI[wrongPos] <- -errorScore
    EnegI[rightPosI] <- -fitScore
    EnegI[rightPosII] <- -fitScore*fitMult
    
  }

  return(list(exprs = NEMlist$exprs, fc = CompMatCont, E0 = E0, Epos = Epos, Eneg = Eneg, EposI = EposI, EnegI = EnegI))
}
