simAnn <- function(CNOlist, NEMlist, model, sapars=c(T,HT,f1,f2,k,Tstop), approach="fc", parameters=c(0.5,1), verbose=FALSE) {

  bitStart <- sample(c(0,1), length(model$reacID), replace = TRUE)

  bitStrings <- bitStart

  bitStringsScores <- 0
}
