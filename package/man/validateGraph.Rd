\name{validateGraph}
\alias{validateGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualize fitted data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Visualizes the data in context of the given hyper-graph and the regulation of each E-genes. A heatmap shows the expected response scheme of the S-gene S and the observed response scheme of the E-genes regulated by S.
}
\usage{
validateGraph(CNOlist, fc=NULL, exprs=NULL, approach = "fc", model, bString, Egenes = 25, Sgene = 1, parameters = list(cutOffs = c(0, 1, 0), scoring = c(0.1, 0.2, 0.9)), plot = TRUE, disc = 0, affyIds = TRUE, sim = 0, relFit = FALSE, complete = FALSE, xrot = 25, Rowv = F, Colv = F, dendrogram = "none", soft = TRUE, colSideColors = NULL, affychip = "hgu133plus2", method = "s", ranks = F, breaks = NULL, col = "RdYlGn", csc = TRUE, sizeFac = 10^-10, verbose = T, order = "rank", colnames = "bio", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
CNOlist object.
}
  \item{fc}{
%%     ~~Describe \code{NEMlist} here~~
matrix of foldchanges (observed response scheme, ORS).
}
  \item{exprs}{
%%     ~~Describe \code{NEMlist} here~~
optional matrix of normalized expression (observed activation scheme).
}
  \item{approach}{
%%     ~~Describe \code{approach} here~~
not used
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Model object.
}
  \item{bString}{
%%     ~~Describe \code{bString} here~~
Binary string denoting the hyper-graph.
}
  \item{Egenes}{
%%     ~~Describe \code{Egenes} here~~
Maximal number of visualized E-genes.
}
  \item{Sgene}{
%%     ~~Describe \code{Sgene} here~~
Integer denoting the S-gene. See colnames(CNOlist@signals[[1]]) to match integer with S-gene name.
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
not used
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
Plot the heatmap. If FALSE, only corresponding information is outputed.
}
  \item{disc}{
%%     ~~Describe \code{disc} here~~
Discretize the data.
}
  \item{affyIds}{
%%     ~~Describe \code{affyIds} here~~
Experimental. Turn Affymetrix Ids into HGNC gene symbols.
}
  \item{sim}{
%%     ~~Describe \code{sim} here~~
not used
}
  \item{relFit}{
%%     ~~Describe \code{relFit} here~~
not used
}
  \item{complete}{
%%     ~~Describe \code{complete} here~~
not used
}
  \item{xrot}{
%%     ~~Describe \code{xrot} here~~
See ?heatmapOP.
}
  \item{Rowv}{
%%     ~~Describe \code{Rowv} here~~
See ?heatmapOP.
}
  \item{Colv}{
%%     ~~Describe \code{Colv} here~~
See ?heatmapOP.
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
See ?heatmapOP.
}
  \item{soft}{
%%     ~~Describe \code{soft} here~~
not used
}
  \item{colSideColors}{
%%     ~~Describe \code{colSideColors} here~~
See ?heatmapOP.
}
  \item{affychip}{
%%     ~~Describe \code{affychip} here~~
Define Affymetrix chip used to generate the data.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
Scoring method can be a correlation or distance measure. See ?cor and ?dist for details.
}
  \item{ranks}{
%%     ~~Describe \code{ranks} here~~
Turn data into ranks.
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
See ?heatmapOP.
}
  \item{col}{
%%     ~~Describe \code{col} here~~
See ?heatmapOP.
}
  \item{csc}{
%%     ~~Describe \code{csc} here~~
not used
}
  \item{sizeFac}{
%%     ~~Describe \code{sizeFac} here~~
Size factor penelizing the hyper-graph size.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Verbose output.
}
  \item{order}{
%%     ~~Describe \code{order} here~~
Order by "rank" or "name".
}
  \item{colnames}{
%%     ~~Describe \code{colnames} here~~
not used
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
\item{v}{
Includes different information like how well each E-gene fits.
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
CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2, signal = c("A", "B","C","D"))
model <- preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
fc <- computeFc(CNOlist, exprs)
initBstring <- rep(0, length(model$reacID))
res <- bnem(search = "greedy", CNOlist = CNOlist, NEMlist = NEMlist, model = model, parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE, maxSteps = Inf)
rownames(fc) <- 1:nrow(fc)
val <- validateGraph(CNOlist = CNOlist, NEMlist = NEMlist, model = model, bString = res$bString, Egenes = 10, Sgene = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
\keyword{ visualization }% __ONLY ONE__ keyword per line
