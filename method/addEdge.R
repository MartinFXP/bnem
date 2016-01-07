addEdge <- function(edges, CNOlist, model, n = 100, full = FALSE) {
  if (full) {
    sifMatrix <- numeric()
    graph <- model$reacID[-grep("\\+", model$reacID)]
    for (i in c(graph, edges)) {
      tmp2 <- unlist(strsplit(i, "="))
      if (gsub("!", "", tmp2[1]) %in% tmp2[1]) {
        sifMatrix <- rbind(sifMatrix, c(tmp2[1], 1, tmp2[2]))
      } else {
        sifMatrix <- rbind(sifMatrix, c(gsub("!", "", tmp2[1]), -1, tmp2[2]))
      }
    }
    write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    PKN2 <- readSIF("temp.sif")
    unlink("temp.sif")
    checkSignals(CNOlist,PKN2)
    verbose <- FALSE
    indices<-indexFinder(CNOlist,PKN2,verbose=verbose)
    NCNOindices<-findNONC(PKN2,indices,verbose=verbose)
    NCNOcut<-cutNONC(PKN2,NCNOindices)
    indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
    NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
    indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
    NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate = n)
    resECNOlist<-residualError(CNOlist)
    Fields4Sim<-prep4sim(NCNOcutCompExp)
    model2 <- NCNOcutCompExp
  } else {
    sifMatrix <- numeric()
    graph <- model$reacID
    for (i in c(graph, edges)) {
      input <- unlist(strsplit(i, "="))
      output <- input[2]
      inputs <- unlist(strsplit(input[1], "\\+"))
      for (j in inputs) {
        if (gsub("!", "", j) %in% j) {
          sifMatrix <- rbind(sifMatrix, c(j, 1, output))
        } else {
          sifMatrix <- rbind(sifMatrix, c(gsub("!", "", j), -1, output))
        }
      }
    }
    write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    PKN2 <- readSIF("temp.sif")
    unlink("temp.sif")
    checkSignals(CNOlist,PKN2)
    indices<-indexFinder(CNOlist,PKN2,verbose=FALSE)
    NCNOindices<-findNONC(PKN2,indices,verbose=FALSE)
    NCNOcut<-cutNONC(PKN2,NCNOindices)
    indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
    NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
    indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
    NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate = n)
    resECNOlist<-residualError(CNOlist)
    Fields4Sim<-prep4sim(NCNOcutCompExp)
    model2 <- NCNOcutCompExp
  }
  return(model2)
}
