\name{simulateStatesRecursive}
\alias{simulateStatesRecursive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate states for a model and CNO object
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simulate the steady state of all vertices of a hyper-graph given a list of stimulated and inhibited vertices.
}
\usage{
simulateStatesRecursive(CNOlist, model, bString, NEMlist = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Model object.
}
  \item{bString}{
%%     ~~Describe \code{bString} here~~
Binary vector.
}
  \item{NEMlist}{
%%     ~~Describe \code{NEMlist} here~~
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
\item{states}{
Matrix including the steady states of all vertices (rows) for all conditions (stimulated and inhibited vertices).
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
model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
states <- simulateStatesRecursive(CNOlist, model, rep(1, length(model$reacID)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
