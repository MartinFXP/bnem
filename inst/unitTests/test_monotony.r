sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
c("C", 1, "D"))
file <- tempfile("interaction",fileext=".sif")
write.table(sifMatrix, file = file, sep = "\t",
row.names = FALSE, col.names = FALSE,
quote = FALSE)
PKN <- CellNOptR::readSIF(file)
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1,
maxInhibit = 2, signals = c("A", "B","C","D"))
model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
nrow(slot(CNOlist, "cues")))
fc <- computeFc(CNOlist, exprs)
initBstring <- rep(0, length(model$reacID))
res <- bnem(search = "greedy", model = model, CNOlist = CNOlist,
fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
maxSteps = Inf)
checkTrue(all(res$scores[[1]][-1]-res$scores[[1]][-length(res$scores[[1]])]<=0))
