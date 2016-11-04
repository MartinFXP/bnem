\name{plotDnf}
\alias{plotDnf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
draw dnf
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDnf(dnf = NULL, freq = NULL, stimuli = c(), signals = c(), inhibitors = c(), connected = TRUE, CNOlist = NULL, cex = NULL, fontsize = NULL, labelsize = NULL, type = 2, lwd = 2, edgelwd = 2, legend = 0, x = 0, y = 0, xjust = 0, yjust = 0, width = 1.5, height = 1, rankdir = "TB", rank = "same", layout = "dot", main = "", sub = "", cex.main = 1.5, cex.sub = 1, col.sub = "grey", fontcolor = NULL, nodestates = NULL, simulate = NULL, andcolor = "transparent", edgecol = NULL, labels = NULL, labelcol = "blue", nodelabel = NULL, nodecol = NULL, bordercol = NULL, nodeshape = NULL, verbose = FALSE, edgestyle = NULL, nodeheight = NULL, nodewidth = NULL, edgewidth = NULL, lty = NULL, hierarchy = NULL, showall = FALSE, nodefontsize = NULL, edgehead = NULL, edgelabel = NULL, edgetail = NULL, bool = TRUE, draw = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnf}{
%%     ~~Describe \code{dnf} here~~
Hyper-graph in disjunctive normal form.
}
  \item{freq}{
%%     ~~Describe \code{freq} here~~
Frequency of hyper-edges.
}
  \item{stimuli}{
%%     ~~Describe \code{stimuli} here~~
Vertices which can be stimulated.
}
  \item{signals}{
%%     ~~Describe \code{signals} here~~
Vertices which regulate E-genes.
}
  \item{inhibitors}{
%%     ~~Describe \code{inhibitors} here~~
Vertices which can be inhibited.
}
  \item{connected}{
%%     ~~Describe \code{connected} here~~
If TRUE, only includes vertices which are connected to other vertices.
}
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object. Optional instead of stimuli, inhibitors or signals.
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
Global font size.
}
  \item{fontsize}{
%%     ~~Describe \code{fontsize} here~~
Vertice label size.
}
  \item{labelsize}{
%%     ~~Describe \code{labelsize} here~~
Edge label size.
}
  \item{type}{
%%     ~~Describe \code{type} here~~
Different plot types.
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
Line width.
}
  \item{edgelwd}{
%%     ~~Describe \code{edgelwd} here~~
Edgeline width.
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
0 shows no legend. 1 shows legend as a graph. 2 shows legend in a box.
}
  \item{x}{
%%     ~~Describe \code{x} here~~
x coordinate of box legend.
}
  \item{y}{
%%     ~~Describe \code{y} here~~
y coordinate of box legend.
}
  \item{xjust}{
%%     ~~Describe \code{xjust} here~~
Justification of legend box left, right or center (-1,1,0).
}
  \item{yjust}{
%%     ~~Describe \code{yjust} here~~
Justification of legend box top, bottom or middle (-1,1,0).
}
  \item{width}{
%%     ~~Describe \code{width} here~~
Vertice width.
}
  \item{height}{
%%     ~~Describe \code{height} here~~
Vertice height.
}
  \item{rankdir}{
%%     ~~Describe \code{rankdir} here~~
not used
}
  \item{rank}{
%%     ~~Describe \code{rank} here~~
not used
}
  \item{layout}{
%%     ~~Describe \code{layout} here~~
Graph layout. See graphvizCapabilities()$layoutTypes.
}
  \item{main}{
%%     ~~Describe \code{main} here~~
Main title.
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
Subtitle.
}
  \item{cex.main}{
%%     ~~Describe \code{cex.main} here~~
Main title font size.
}
  \item{cex.sub}{
%%     ~~Describe \code{cex.sub} here~~
Subtitle font size.
}
  \item{col.sub}{
%%     ~~Describe \code{col.sub} here~~
Font color of subtitle.
}
  \item{fontcolor}{
%%     ~~Describe \code{fontcolor} here~~
Global font color.
}
  \item{nodestates}{
%%     ~~Describe \code{nodestates} here~~
Binary state of each vertice.
}
  \item{simulate}{
%%     ~~Describe \code{simulate} here~~
Simulate stimulation and inhibition of a list of vertices. E.g. simulate = list(stimuli = c("A", "B"), inhibitors = c("C", "D")).
}
  \item{andcolor}{
%%     ~~Describe \code{andcolor} here~~
not used
}
  \item{edgecol}{
%%     ~~Describe \code{edgecol} here~~
Vector with colors for every edge of the graph (not hyper-graph). E.g. an AND gate consists of three distinct edges.
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
Vector with labels for the edges.
}
  \item{labelcol}{
%%     ~~Describe \code{labelcol} here~~
Vector with label colors for the edges.
}
  \item{nodelabel}{
%%     ~~Describe \code{nodelabel} here~~
List of vertices with labels as input. E.g. labels = list(A="test", B="label for B").
}
  \item{nodecol}{
%%     ~~Describe \code{nodecol} here~~
List of vertices with colors as input.
}
  \item{bordercol}{
%%     ~~Describe \code{bordercol} here~~
List of vertices with colors as input.
}
  \item{nodeshape}{
%%     ~~Describe \code{nodeshape} here~~
List of vertices with shapes (diamond, box, square,...).
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Verbose output.
}
  \item{edgestyle}{
%%     ~~Describe \code{edgestyle} here~~
not used
}
  \item{nodeheight}{
%%     ~~Describe \code{nodeheight} here~~
List of vertices with height as input.
}
  \item{nodewidth}{
%%     ~~Describe \code{nodewidth} here~~
List of vertices with width as input.
}
  \item{edgewidth}{
%%     ~~Describe \code{edgewidth} here~~
Vector with edge widths.
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
Vector with edge styles (line, dotted,...).
}
  \item{hierarchy}{
%%     ~~Describe \code{hierarchy} here~~
List with the hierarchy of the vertices. E.g. list(top = c("A", "B"), bottom = c("C", "D")).
}
  \item{showall}{
%%     ~~Describe \code{showall} here~~
See "connected" above.
}
  \item{nodefontsize}{
%%     ~~Describe \code{nodefontsize} here~~
not used
}
  \item{edgehead}{
%%     ~~Describe \code{edgehead} here~~
Vector with edge heads.
}
  \item{edgelabel}{
%%     ~~Describe \code{edgelabel} here~~
Vector with edge labels.
}
  \item{edgetail}{
%%     ~~Describe \code{edgetail} here~~
Vector with edge tails.
}
  \item{bool}{
%%     ~~Describe \code{bool} here~~
If TRUE, only shows normal graph and no AND gates.
}
  \item{draw}{
%%     ~~Describe \code{draw} here~~
Do not plot the graph and only output the graphNEL object.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
g <- c("A=B", "C+D=E", "!G=H")
g2 <- plotDnf(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
