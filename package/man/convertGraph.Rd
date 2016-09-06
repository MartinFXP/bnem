\name{convertGraph}
\alias{convertGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
distributivity for disjunctive normal form (dnf) to conjunctive normal form (cnf) or vice versa
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Turns a dnf into a cnf and vice versa by the usage of the distributivity laws.
}
\usage{
convertGraph(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
Dnf or cnf.
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
\item{h}{
%%     ~~Describe \code{g} here~~
Dnf or cnf.
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
g <- c("A+B=E", "C+D=E")
h <- convertGraph(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
