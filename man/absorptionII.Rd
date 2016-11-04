\name{absorptionII}
\alias{absorptionII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
applies a third absoprtion law to a dnf
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A Boolean hyper-graph is reduced by a third absorption law.
}
\usage{
absorptionII(bString, model = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bString}{
%%     ~~Describe \code{bString} here~~
A binary vector denoting which hyper-edges are present or a graph in disjunctive normal form.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
If bString is a binary vector, a model including the complete Boolean hyper-graph space.
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
A reduced binary vector or hyper-graph.
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
\code{\link{absorption}}
}
\examples{
graph <- c("A+B=C", "A=C")
absorptionII(graph)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Boolean algebra }
