\name{transRed}
\alias{transRed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Transitive reduction of a Boolean graph.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the transitive reduction of a hyper-graph in disjunctive normal form.
}
\usage{
transRed(g, max.iter = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
Hyper-graph in dnf.
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
Transitive reduction of g.
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
g <- c("A=B", "A=C", "B=C", "B=D", "!A=D")
gred <- transRed(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
