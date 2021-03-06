% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bnem_main.r
\name{plot.bnemBs}
\alias{plot.bnemBs}
\title{Plot Bootstrap result}
\usage{
\method{plot}{bnemBs}(
  x,
  scale = 3,
  shift = 0.1,
  cut = 0.5,
  dec = 2,
  ci = 0,
  cip = 0.95,
  method = "exact",
  ...
)
}
\arguments{
\item{x}{bnemBs object}

\item{scale}{numeric value for scaling the edgewidth}

\item{shift}{numeric value for shifting the edgewidth}

\item{cut}{shows only edges with a fraction larger than cut}

\item{dec}{integer for function round}

\item{ci}{if TRUE shows confidence intervals}

\item{cip}{range for the confidence interval, e.g. 0.95}

\item{method}{method to use for conidence interval
computation (see function binom.confint from package binom)}

\item{...}{additional parameters for the function mnem::plotDnf}
}
\value{
plot of the network from the bootstrap
}
\description{
Shows the result of a Boostrap with either edge frequencies
or confidence intervals
}
\examples{
sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
c("C", 1, "D"))
temp.file <- tempfile(pattern="interaction",fileext=".sif")
write.table(sifMatrix, file = temp.file, sep = "\t",
row.names = FALSE, col.names = FALSE,
quote = FALSE)
PKN <- CellNOptR::readSIF(temp.file)
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1,
maxInhibit = 2, signals = c("A", "B","C","D"))
model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
expression <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
nrow(slot(CNOlist, "cues")))
fc <- computeFc(CNOlist, expression)
initBstring <- rep(0, length(model$reacID))
res <- bnemBs(search = "greedy", model = model, CNOlist = CNOlist,
fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
maxSteps = Inf)
}
\author{
Martin Pirkl
}
