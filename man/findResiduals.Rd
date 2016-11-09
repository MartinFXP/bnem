\name{findResiduals}
\alias{findResiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate and visualize local residuals of the network compared to data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a network hypothesis and data, this function calculates residuals
caused by the ERS of the network.
}
\usage{
findResiduals(bString, CNOlist, model, fc=NULL, exprs=NULL, egenes=NULL, NEMlist = NULL, parameters = list(cutOffs = c(0,1,0), scoring = c(0.1,0.2,0.9)), method = "s", sizeFac = 10^-10, main = "residuals for decoupled vertices", sub = "green residuals are added effects (left positive, right negative) and red residuals are deleted effects", cut = TRUE, approach = "fc", parallel = NULL, verbose = TRUE)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{bString}{
%%     ~~Describe \code{vertices} here~~
Binary vector denoting the network given a model.
}
  \item{CNOlist}{
%%     ~~Describe \code{negation} here~~
CNOlist object.
}
  \item{model}{
%%     ~~Describe \code{max.edge.size} here~~
Network model object.
}
  \item{fc}{
%%     ~~Describe \code{max.edge.size} here~~
ORS of the data.
}
  \item{exprs}{
%%     ~~Describe \code{max.edge.size} here~~
Optional activation scheme of the data.
}
  \item{egenes}{
%%     ~~Describe \code{max.edge.size} here~~
Atachment of the E-genes (optional).
}
  \item{NEMlist}{
%%     ~~Describe \code{max.edges} here~~
NEMlist object (optional).
}
  \item{parameters}{
%%     ~~Describe \code{max.edge.size} here~~
not used
}
  \item{method}{
%%     ~~Describe \code{max.edge.size} here~~
Scoring method (optional).
}
  \item{sizeFac}{
%%     ~~Describe \code{max.edge.size} here~~
Zeta parameter to penelize netowkr size.
}
  \item{main}{
%%     ~~Describe \code{max.edge.size} here~~
Main title of the figure.
}
  \item{sub}{
%%     ~~Describe \code{max.edge.size} here~~
Subtitle of the figure.
}
  \item{cut}{
%%     ~~Describe \code{max.edge.size} here~~
If TRUE does not visualize experiments/S-genes which do nto have any residuals.
}
  \item{approach}{
%%     ~~Describe \code{max.edge.size} here~~
not used
}
  \item{parallel}{
%%     ~~Describe \code{max.edge.size} here~~
Set number of cores/threads.
}
  \item{verbose}{
%%     ~~Describe \code{max.edge.size} here~~
verbose output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{residuals}{
Matrices and figure denoting potential residuals in the network.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Martin Pirkl
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(bnem)
     
     sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
     write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
     quote = FALSE)
     PKN <- readSIF("temp.sif")
     unlink('temp.sif')
     CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
     model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
     exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
     fc <- computeFc(CNOlist, exprs)
     initBstring <- rep(0, length(model$reacID))
     res <- bnem(search = "greedy", CNOlist = CNOlist, fc = fc, model = model, parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
     rownames(fc) <- 1:nrow(fc)
     ## val <- validateGraph(CNOlist = CNOlist, fc = fc, model = model, bString = res$bString, Egenes = 10, Sgene = 4)
     residuals <- findResiduals(res$bString, CNOlist, model, fc = fc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
