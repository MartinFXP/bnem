\name{transClose}
\alias{transClose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate the transitive closure of a graph
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes the transitive closure of a Boolean hyper-graph.
}
\usage{
transClose(g, max.iter = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
Hyper-graph in its disjunctive normal form.
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
Maximal iterations. If the algorithm does not converge, the output is only an approximation.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Verbose output.
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
Hyper-graph in dnf depicting the transitive closure of g.
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
g <- c("A=B", "B=C")
gclose <- transClose(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
