\name{bnem-package}
\alias{bnem-package}
\alias{bnem}
\docType{package}
\title{
Boolean Nested Effects Models.
}
\description{
Estimates a signalling pathway from gene expression data derived from perturbation experiments.
}
\details{
\tabular{ll}{
Package: \tab bnem\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-04-28\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Martin Pirkl

Maintainer: Martin Pirkl <martin-franz-xaver.pirkl@ukr.de>
~~ The author and/or maintainer of the package ~~
}
\references{
Pirkl, Martin, Hand, Elisabeth, Kube, Dieter, & Spang, Rainer. 2016. Analyzing synergistic and non-synergistic interactions in signalling pathways using Boolean Nested Effect Models. Bioinformatics, 32(6), 893–900.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(bnem)

sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
quote = FALSE)
PKN <- readSIF("temp.sif")
unlink('temp.sif')
model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
NEMlist <- list()
NEMlist$exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
NEMlist$fc <- computeFc(CNOlist, NEMlist$exprs)
initBstring <- rep(0, length(model$reacID))
res <- bnem(search = "greedy", CNOlist = CNOlist, NEMlist = NEMlist, model = model, parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
}
