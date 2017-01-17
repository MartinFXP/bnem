#' @noRd
preprocessInput <- function(stimuli=NULL, inhibitors=NULL, signals=NULL, design = NULL, exprs=NULL, fc=NULL, pkn, maxInputsPerGate=100) {

    if (is.null(design)) {

        stimcols <- stimcols2 <- matrix(0, length(stimuli), ncol(fc))

        for (i in 1:length(stimuli)) {

            tmp <- numeric(ncol(fc))

            tmp[grep(stimuli[i], gsub("_vs_.*", "", colnames(fc)))] <- 1

            stimcols[i, ] <- tmp

            tmp <- numeric(ncol(fc))

            tmp[grep(stimuli[i], gsub(".*_vs_", "", colnames(fc)))] <- 1

            stimcols2[i, ] <- tmp

        }

        maxStim <- max(c(apply(stimcols, 2, sum), apply(stimcols2, 2, sum)))

        inhibitcols <- inhibitcols2 <- matrix(0, length(inhibitors), ncol(fc))

        for (i in 1:length(inhibitors)) {

            tmp <- numeric(ncol(fc))

            tmp[grep(inhibitors[i], gsub("_vs_.*", "", colnames(fc)))] <- 1

            inhibitcols[i, ] <- tmp

            tmp <- numeric(ncol(fc))

            tmp[grep(inhibitors[i], gsub(".*_vs_", "", colnames(fc)))] <- 1

            inhibitcols2[i, ] <- tmp

        }

        maxInhibit <- max(c(apply(inhibitcols, 2, sum), apply(inhibitcols2, 2, sum)))

        if (is.null(signals)) { signals <- unique(c(stimuli, inhibitors)) }

        CNOlist <- dummyCNOlist(stimuli=stimuli, inhibitors=inhibitors, maxStim=maxStim, maxInhibit=maxInhibit, signals=signals)

        model <- preprocessing(CNOlist, pkn, maxInputsPerGate=maxInputsPerGate)

    }

    NEMlist <- list()

    NEMlist$fc <- fc

    if (is.null(exprs)) {
        NEMlist$exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10, nrow(CNOlist@cues))
    } else {
        NEMlist$exprs <- exprs
    }

    return(list(CNOlist=CNOlist, model=model, NEMlist=NEMlist))

}



