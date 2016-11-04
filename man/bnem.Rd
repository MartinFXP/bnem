\name{bnem}
\alias{bnem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
main function to learn the Boolean network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The main high level function to learn an optimal Boolean hyper-graph describing the data.
}
\usage{
bnem(search = "loc", fc=NULL, exprs=NULL, pkn, design=NULL,
stimuli=NULL, inhibitors=NULL, signals=NULL, CNOlist=MULL, model=NULL, sizeFac = 10^-10, NAFac = 1, NEMlist, parameters = list(cutOffs = c(0, 1, 0), scoring = c(0.1, 0.2, 0.9)), parallel = NULL, method = "s", relFit = F, verbose = TRUE, reduce = TRUE, initBstring = NULL, popSize = 100, pMutation = 0.5, maxTime = Inf, maxGens = Inf, stallGenMax = 10, relTol = 0.01, priorBitString = NULL, selPress = c(1.2, 1e-04), approach = "fc", fit = "linear", targetBstring = "none", elitism = NULL, inversion = NULL, parallel2 = 1, selection = c("t"), type = "SOCK", exhaustive = FALSE, delcyc = TRUE, seeds = 1, maxSteps = Inf, node = NULL, absorpII = TRUE, draw = TRUE, prior = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{search}{
%%     ~~Describe \code{search} here~~
Type of search heuristic. Either "greedy", "genetic" or "exhaustive". "greedy" uses a greedy algorithm to move through the local neighbourhood of a initial hyper-graph. "genetic" uses a genetic algorithm. "exhaustive" searches through the complete search space and is not recommended.
}
  \item{fc}{
%%     ~~Describe \code{CNOlist} here~~
Foldchanges of gene expression values or equivalent input (normalized
pvalues, logodds, ...). If left NULL, the gene expression data is used
to calculate naive foldchanges.
}
  \item{exprs}{
%%     ~~Describe \code{CNOlist} here~~
Optional normalized gene expression data.
}
  \item{pkn}{
%%     ~~Describe \code{CNOlist} here~~
Prior knowledge network.
}
  \item{design}{
%%     ~~Describe \code{CNOlist} here~~
Optional design matrix for the gene expression values, if available. If
kept NULL, bnem need either stimuli, inhibitors or a CNOlist object.
}
  \item{stimuli}{
%%     ~~Describe \code{CNOlist} here~~
Character vector of stimuli names.
}
  \item{inhibitors}{
%%     ~~Describe \code{CNOlist} here~~
Character vector of inhibitors.
}
  \item{signals}{
%%     ~~Describe \code{CNOlist} here~~
Optional character vector of signals. Signals are S-genes, which can
directly regulate E-genes. If left NULL, alls stimuli and inhibitors are
defined as signals.
}
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object if available.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Model object including the search space, if available.
}
  \item{sizeFac}{
%%     ~~Describe \code{sizeFac} here~~
Size factor penelizing the hyper-graph size.
}
  \item{NAFac}{
%%     ~~Describe \code{NAFac} here~~
NE factor penelizing NAs in the data.
}
  \item{NEMlist}{
%%     ~~Describe \code{NEMlist} here~~
NEMlist object including the contrasts.
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
not used
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
Parallelize the search. An integer value specifies the number of threads on the local machine. A list object as in list(c(1,2,3), c("machine1", "machine2", "machine3")) specifies the threads distributed on different machines (local or others).
}
  \item{method}{
%%     ~~Describe \code{method} here~~
Scoring method can be a correlation or distance measure. See ?cor and ?dist for details.
}
  \item{relFit}{
%%     ~~Describe \code{relFit} here~~
not used
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
TRUE gives additional information during the search.
}
  \item{reduce}{
%%     ~~Describe \code{reduce} here~~
not used
}
  \item{initBstring}{
%%     ~~Describe \code{initBstring} here~~
Binary string of the initial hyper-graph.
}
  \item{popSize}{
%%     ~~Describe \code{popSize} here~~
Popultaion size (only "genetic").
}
  \item{pMutation}{
%%     ~~Describe \code{pMutation} here~~
Probability for mutation (only "genetic").
}
  \item{maxTime}{
%%     ~~Describe \code{maxTime} here~~
Define a maximal time for the search.
}
  \item{maxGens}{
%%     ~~Describe \code{maxGens} here~~
Maximal number of generations (only "genetic").
}
  \item{stallGenMax}{
%%     ~~Describe \code{stallGenMax} here~~
Maximum number of stall generations (only "genetic").
}
  \item{relTol}{
%%     ~~Describe \code{relTol} here~~
Score tolerance for networks defined as optimal but with a lower score as the real optimum (only "genetic").
}
  \item{priorBitString}{
%%     ~~Describe \code{priorBitString} here~~
Binary string defining hyper-edges which are added to every hyper-graph. E.g. if you know hyper-edge 55 is definitly there and want to fix that set priorBitString[55] <- 1 (only "genetic").
}
  \item{selPress}{
%%     ~~Describe \code{selPress} here~~
Selection pressure for the stochastic universal sampling (only "genetic").
}
  \item{approach}{
%%     ~~Describe \code{approach} here~~
not used
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
not used
}
  \item{targetBstring}{
%%     ~~Describe \code{targetBstring} here~~
not used
}
  \item{elitism}{
%%     ~~Describe \code{elitism} here~~
Number of best hyper-graph transferred to the next generation  (only "genetic").
}
  \item{inversion}{
%%     ~~Describe \code{inversion} here~~
Number of worst hyper-graphs for which their binary strings are inversed  (only "genetic").
}
  \item{parallel2}{
%%     ~~Describe \code{parallel2} here~~
not used
}
  \item{selection}{
%%     ~~Describe \code{selection} here~~
"t" for tournament selection and "s" for stochastic universal sampling (only "genetic").
}
  \item{type}{
%%     ~~Describe \code{type} here~~
not used
}
  \item{exhaustive}{
%%     ~~Describe \code{exhaustive} here~~
If TRUE an exhaustive search is conducted if the genetic algorithm would take longer (only "genetic").
}
  \item{delcyc}{
%%     ~~Describe \code{delcyc} here~~
If TRUE deleted cycles in all hyper-graphs (recommended).
}
  \item{seeds}{
%%     ~~Describe \code{seeds} here~~
not used
}
  \item{maxSteps}{
%%     ~~Describe \code{maxSteps} here~~
Maximal number of steps (only "greedy").
}
  \item{node}{
%%     ~~Describe \code{node} here~~
not used
}
  \item{absorpII}{
%%     ~~Describe \code{absorpII} here~~
Use a thrid absorption law.
}
  \item{draw}{
%%     ~~Describe \code{draw} here~~
If TRUE draws the network evolution.
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
Binary vector. A 1 specifies hyper-edges which should not be optimized (only "greedy").
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
not used
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
\item{res}{
List object including the optimized hyper-graph and its corresponding binary string.
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
sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"), c("C", 1, "D"))
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE, col.names = FALSE,
quote = FALSE)
PKN <- readSIF("temp.sif")
unlink('temp.sif')
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signals = c("A", "B","C","D"))
model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
fc <- computeFc(CNOlist, exprs)
initBstring <- rep(0, length(model$reacID))
res <- bnem(search = "greedy", model = model, CNOlist = CNOlist, fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"), parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathway }
\keyword{ optimization }
\keyword{ network }% __ONLY ONE__ keyword per line
