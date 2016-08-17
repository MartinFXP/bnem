\name{heatmapOP}
\alias{heatmapOP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
draw heatmap
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses the functionality of the lattice package to draw a heatmap of a matrix.
}
\usage{
heatmapOP(x, col = "RdYlGn", coln = 11, bordercol = "grey", borderwidth = 0.1, breaks = "sym", main = "heatmap by Oscar Perpi\303\261\303\241n", sub = "http://oscarperpinan.github.io/rastervis/; http://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps", dendrogram = "none", colorkey = list(space = "right"), Colv = TRUE, Rowv = TRUE, xrot = 90, yrot = 0, shrink = c(1, 1), cexCol = 1, cexRow = 1, cexMain = 1, cexSub = 1, colSideColors = NULL, aspect = "fill", contour = FALSE, useRaster = FALSE, xlab = NULL, ylab = NULL, colSideColorsPos = "top", clust = NULL, clusterx = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Matrix.
}
  \item{col}{
%%     ~~Describe \code{col} here~~
Color. See brewer.pal.info for all available color schemes.
}
  \item{coln}{
%%     ~~Describe \code{coln} here~~
Number of colors.
}
  \item{bordercol}{
%%     ~~Describe \code{bordercol} here~~
Border color.
}
  \item{borderwidth}{
%%     ~~Describe \code{borderwidth} here~~
Border width.
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
Defines the breaks in the color range. "sym" makes the breaks symmetric around 0.
}
  \item{main}{
%%     ~~Describe \code{main} here~~
Main title.
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
Subtitle.
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
Draw dendrogram with "both", "col" or "row", or do not draw with "none".
}
  \item{colorkey}{
%%     ~~Describe \code{colorkey} here~~
Draw colorkey list(space="left") or list(space="right"). 
}
  \item{Colv}{
%%     ~~Describe \code{Colv} here~~
Cluster columns (TRUE) or not (FALSE).
}
  \item{Rowv}{
%%     ~~Describe \code{Rowv} here~~
Cluster rows (TRUE) or not (FALSE).
}
  \item{xrot}{
%%     ~~Describe \code{xrot} here~~
Rotate the column names by degree.
}
  \item{yrot}{
%%     ~~Describe \code{yrot} here~~
Rotate the row names by degree.
}
  \item{shrink}{
%%     ~~Describe \code{shrink} here~~
c(x,y) defines a range of size for the data boxes from low to high.
}
  \item{cexCol}{
%%     ~~Describe \code{cexCol} here~~
Font size of column names.
}
  \item{cexRow}{
%%     ~~Describe \code{cexRow} here~~
Font size of row names.
}
  \item{cexMain}{
%%     ~~Describe \code{cexMain} here~~
Font size of main title.
}
  \item{cexSub}{
%%     ~~Describe \code{cexSub} here~~
Font size of subtitle.
}
  \item{colSideColors}{
%%     ~~Describe \code{colSideColors} here~~
Defines a numeric vector to annotate columns with different colors.
}
  \item{aspect}{
%%     ~~Describe \code{aspect} here~~
"iso" for quadratic boxes or "fill" for streched boxes.
}
  \item{contour}{
%%     ~~Describe \code{contour} here~~
TRUE adds a contour plot.
}
  \item{useRaster}{
%%     ~~Describe \code{useRaster} here~~
not used
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
Label for the x-axis.
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
Label for the y-axis.
}
  \item{colSideColorsPos}{
%%     ~~Describe \code{colSideColorsPos} here~~
Place colSideColors at the "top" or "bottom".
}
  \item{clust}{
%%     ~~Describe \code{clust} here~~
not used
}
  \item{clusterx}{
%%     ~~Describe \code{clusterx} here~~
Optional data matrix y with the same dimensions as x. x is columns or rows are sorted by the cluster information of y.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Optional arguments.
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Martin Pirkl and Oscar Perpiñán @ http://oscarperpinan.github.io/rastervis/
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- matrix(rnorm(50), 10, 5)
heatmapOP(x, dendrogram = "both", aspect = "iso", xrot = 45)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
