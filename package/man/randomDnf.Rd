\name{randomDnf}
\alias{randomDnf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create random dnf
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Create a random hyper-graph in its disjunctive normal form with different restrictions and features.
}
\usage{
randomDnf(vertices = 10, negation = TRUE, max.edge.size = NULL, max.edges = NULL, dag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vertices}{
%%     ~~Describe \code{vertices} here~~
Maximal number of vertices.
}
  \item{negation}{
%%     ~~Describe \code{negation} here~~
TRUE allows negation.
}
  \item{max.edge.size}{
%%     ~~Describe \code{max.edge.size} here~~
Maximal number of parents per hyper-edge.
}
  \item{max.edges}{
%%     ~~Describe \code{max.edges} here~~
Maximal number of hyper-edges.
}
  \item{dag}{
%%     ~~Describe \code{dag} here~~
TRUE results in a acyclic graph.
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
\item{g}{
Dnf describing a hyper-graph.
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
g <- randomDnf(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
