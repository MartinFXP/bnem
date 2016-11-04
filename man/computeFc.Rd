\name{computeFc}
\alias{computeFc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute foldchanges from expression data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in a matrix with specific annotation in the colnames and returns pairwise contrasts between different columns.
}
\usage{
computeFc(CNOlist, y, test = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object.
}
  \item{y}{
%%     ~~Describe \code{y} here~~
Expression matrix.
}
  \item{test}{
%%     ~~Describe \code{test} here~~
For debugging.
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
\item{x}{
Matrix with pairwise contrasts.
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
stimuli <- c("A", "B")
inhibitors <- c("C", "D")
CNOlist <- dummyCNOlist(stimuli, inhibitors, maxStim = 1, maxInhibit = 1, signals = c(stimuli, inhibitors))
exprs <- matrix(rnorm(100*nrow(CNOlist@cues)), 100, nrow(CNOlist@cues))
contrasts <- computeFc(CNOlist, y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
