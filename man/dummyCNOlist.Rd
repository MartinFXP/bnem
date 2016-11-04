\name{dummyCNOlist}
\alias{dummyCNOlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a dummy CNO object
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dummyCNOlist(stimuli = NULL, inhibitors = NULL, maxStim = 0, maxInhibit = 0, signals = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stimuli}{
%%     ~~Describe \code{stimuli} here~~
}
  \item{inhibitors}{
%%     ~~Describe \code{inhibitors} here~~
}
  \item{maxStim}{
%%     ~~Describe \code{maxStim} here~~
}
  \item{maxInhibit}{
%%     ~~Describe \code{maxInhibit} here~~
}
  \item{signals}{
%%     ~~Describe \code{signals} here~~
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (stimuli = NULL, inhibitors = NULL, maxStim = 0, maxInhibit = 0, 
    signals = NULL) 
{
    if (is.null(signals)) {
        signals <- c(stimuli, inhibitors)
    }
    stimn <- length(stimuli)
    inhibn <- length(inhibitors)
    if (maxStim > stimn) {
        maxStim <- stimn
    }
    if (stimn > 0 & maxStim > 0) {
        mat.size <- 0
        for (i in 1:maxStim) {
            mat.size <- mat.size + choose(stimn, i)
        }
        stimDesign <- matrix(0, mat.size, stimn)
        if (length(stimuli) == 1) {
            stimDesign <- t(stimDesign)
        }
        diag(stimDesign) <- 1
        count <- stimn
        if (maxStim > 1) {
            for (i in 2:maxStim) {
                combs <- t(combn(1:stimn, i))
                for (j in 1:nrow(combs)) {
                  count <- count + 1
                  stimDesign[count, combs[j, ]] <- 1
                }
            }
        }
        colnames(stimDesign) <- stimuli
    }
    inhibn <- length(inhibitors)
    if (maxInhibit > inhibn) {
        maxInhibit <- inhibn
    }
    if (inhibn > 0 & maxInhibit > 0) {
        mat.size <- 0
        for (i in 1:maxInhibit) {
            mat.size <- mat.size + choose(inhibn, i)
        }
        inhibDesign <- matrix(0, mat.size, inhibn)
        if (length(inhibitors) == 1) {
            inhibDesign <- t(inhibDesign)
        }
        diag(inhibDesign) <- 1
        count <- inhibn
        if (maxInhibit > 1) {
            for (i in 2:maxInhibit) {
                combs <- t(combn(1:inhibn, i))
                for (j in 1:nrow(combs)) {
                  count <- count + 1
                  inhibDesign[count, combs[j, ]] <- 1
                }
            }
        }
        colnames(inhibDesign) <- inhibitors
    }
    if (stimn > 0 & inhibn > 0) {
        if (maxStim > 0 & maxInhibit > 0) {
            design <- matrix(0, nrow(stimDesign) * nrow(inhibDesign), 
                stimn + inhibn)
            for (i in 1:nrow(stimDesign)) {
                design[((i - 1) * nrow(inhibDesign) + 1):(i * 
                  nrow(inhibDesign)), ] <- cbind(stimDesign[rep(i, 
                  nrow(inhibDesign)), ], inhibDesign)
            }
        }
        if (maxStim > 0 & maxInhibit == 0) {
            inhibDesign <- matrix(0, nrow(stimDesign), inhibn)
            if (length(inhibitors) == 1) {
                inhibDesign <- t(inhibDesign)
            }
            design <- cbind(stimDesign, inhibDesign)
        }
        if (maxStim == 0 & maxInhibit > 0) {
            stimDesign <- matrix(0, nrow(inhibDesign), stimn)
            if (length(stimuli) == 1) {
                stimDesign <- t(stimDesign)
            }
            design <- cbind(stimDesign, inhibDesign)
        }
        if (maxStim == 0 & maxInhibit == 0) {
            inhibDesign <- matrix(0, 1, inhibn)
            stimDesign <- matrix(0, 1, stimn)
            design <- cbind(stimDesign, inhibDesign)
        }
        colnames(design) <- c(stimuli, inhibitors)
    }
    colnamesdesign <- colnames(design)
    design <- rbind(cbind(stimDesign, matrix(0, nrow(stimDesign), 
        (ncol(design) - ncol(stimDesign)))), cbind(matrix(0, 
        nrow(inhibDesign), (ncol(design) - ncol(inhibDesign))), 
        inhibDesign), design)
    colnames(design) <- colnamesdesign
    signaln <- length(signals)
    if (signaln > 0) {
        signalData <- matrix(0, nrow(design) + 1, signaln)
        colnames(signalData) <- signals
    }
    else {
        signalData <- matrix(0, nrow(design) + 1, 1)
    }
    smult <- nrow(design)/nrow(stimDesign)
    imult <- nrow(design)/nrow(inhibDesign)
    design <- rbind(0, design)
    stimDesign <- as.matrix(design[, 1:ncol(stimDesign)])
    inhibDesign <- as.matrix(design[, (ncol(stimDesign) + 1):ncol(design)])
    rownames(design) <- rownames(inhibDesign) <- rownames(stimDesign) <- rownames(signalData) <- c("Ctrl", 
        2:nrow(design))
    getRowname <- function(i, M) {
        r <- paste(colnames(M)[which(M[i, ] == 1)], collapse = "_")
        return(r)
    }
    rownames(design)[2:nrow(design)] <- rownames(inhibDesign)[2:nrow(design)] <- rownames(stimDesign)[2:nrow(design)] <- rownames(signalData)[2:nrow(design)] <- unlist(lapply(as.list(2:nrow(design)), 
        getRowname, design))
    if (ncol(stimDesign) == 1) {
        colnames(stimDesign) <- stimuli
    }
    cnolist <- new("CNOlist", cues = design, inhibitors = inhibDesign, 
        stimuli = stimDesign, signals = list(signalData, signalData), 
        timepoints = as.character(c(0, 1)))
    cnolist <- checkCNOlist(cnolist)
    return(cnolist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
