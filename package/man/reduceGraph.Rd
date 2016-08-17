\name{reduceGraph}
\alias{reduceGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reduce graph
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reduces a given hyper-graph using the absorption law and deleting edges from "dead" vertices.
}
\usage{
reduceGraph(bString, model, CNOlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bString}{
%%     ~~Describe \code{bString} here~~
Binary string denoting the present hyper-edges with 1.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Model object describing the hyper-graph space.
}
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object.
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
\item{bString}{
Reduced binary string.
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
bString <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
