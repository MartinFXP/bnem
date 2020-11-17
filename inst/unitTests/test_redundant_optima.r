set.seed(8675309)
sim <- simBoolGtn()
PKN <- sim$PKN
CNOlist <- sim$CNOlist
model <- sim$model
exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
                nrow(slot(CNOlist, "cues")))
fc <- computeFc(CNOlist, exprs)
initBstring <- rep(0, length(model$reacID))
res <- bnem(search = "genetic", model = model, CNOlist = CNOlist,
            fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
            parallel = NULL, initBstring = initBstring, draw = TRUE,
            verbose = FALSE, maxSteps = Inf)
all.diff <- unlist(lapply(seq_len(nrow(res$bStrings)),function(x) {
    y <- x - res$bString
    y <- all(y==0)
    return(y)
}))
checkTrue(all(!all.diff))
