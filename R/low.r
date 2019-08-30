#' @noRd
preprocessInput <- function(stimuli=NULL, inhibitors=NULL, signals=NULL,
                            design = NULL, exprs=NULL, fc=NULL, pkn,
                            maxInputsPerGate=100) {

    if (is.null(design)) {

        stimcols <- stimcols2 <- matrix(0, length(stimuli), ncol(fc))

        for (i in seq_len(length(stimuli))) {

            tmp <- numeric(ncol(fc))

            tmp[grep(stimuli[i], gsub("_vs_.*", "", colnames(fc)))] <- 1

            stimcols[i, ] <- tmp

            tmp <- numeric(ncol(fc))

            tmp[grep(stimuli[i], gsub(".*_vs_", "", colnames(fc)))] <- 1

            stimcols2[i, ] <- tmp

        }

        maxStim <- max(c(apply(stimcols, 2, sum), apply(stimcols2, 2, sum)))

        inhibitcols <- inhibitcols2 <- matrix(0, length(inhibitors), ncol(fc))

        for (i in seq_len(length(inhibitors))) {

            tmp <- numeric(ncol(fc))

            tmp[grep(inhibitors[i], gsub("_vs_.*", "", colnames(fc)))] <- 1

            inhibitcols[i, ] <- tmp

            tmp <- numeric(ncol(fc))

            tmp[grep(inhibitors[i], gsub(".*_vs_", "", colnames(fc)))] <- 1

            inhibitcols2[i, ] <- tmp

        }

        maxInhibit <- max(c(apply(inhibitcols, 2, sum), apply(inhibitcols2, 2,
                                                              sum)))

        if (is.null(signals)) { signals <- unique(c(stimuli, inhibitors)) }

        CNOlist <- dummyCNOlist(stimuli=stimuli, inhibitors=inhibitors,
                                maxStim=maxStim, maxInhibit=maxInhibit,
                                signals=signals)

        model <- preprocessing(CNOlist, pkn, maxInputsPerGate=maxInputsPerGate)

    }

    NEMlist <- list()

    NEMlist$fc <- fc

    if (is.null(exprs)) {
        NEMlist$exprs <- matrix(rnorm(nrow(CNOlist@cues)*10), 10,
                                nrow(CNOlist@cues))
    } else {
        NEMlist$exprs <- exprs
    }

    return(list(CNOlist=CNOlist, model=model, NEMlist=NEMlist))

}
#' @importFrom graphics abline axis legend mtext par screen split.screen
#' @importFrom utils combn flush.console write.table
#' @noRd
addEdge <-
    function(edges, CNOlist, model, n = 100, full = FALSE) {
        sifMatrix <- numeric()
        graph <- model$reacID[-grep("\\+", model$reacID)]
        for (i in c(graph, edges)) {
            tmp2 <- unlist(strsplit(i, "="))
            if (gsub("!", "", tmp2[1]) %in% tmp2[1]) {
                sifMatrix <- rbind(sifMatrix, c(tmp2[1], 1, tmp2[2]))
            } else {
                sifMatrix <- rbind(sifMatrix, c(gsub("!", "", tmp2[1]), -1,
                                                tmp2[2]))
            }
        }
        write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
        PKN2 <- readSIF("temp.sif")
        unlink("temp.sif")
        model2 <- preprocessing(CNOlist, PKN2, maxInputsPerGate=n)
        if (!full) {
            index <- which(!(model2$reacID %in% model$reacID) &
                           !(model2$reacID %in% edges))
            model2$reacID <- model2$reacID[-index]
            model2$interMat[, -index]
            model2$notMat[, -index]
        }
        return(model2)
    }
#' @noRd
addSignal <-
    function(s, CNOlist, stim = NULL, inhibit = NULL) {
        CNOlist2 <- CNOlist
        CNOlist2@signals[[1]] <- cbind(CNOlist2@signals[[1]], s = 0)
        CNOlist2@signals[[2]] <- cbind(CNOlist2@signals[[2]], s = 0)
        colnames(CNOlist2@signals[[1]]) <- c(colnames(CNOlist@signals[[1]]), s)
        colnames(CNOlist2@signals[[2]]) <- c(colnames(CNOlist@signals[[2]]), s)
        if (!is.null(inhibit)) {
            CNOlist2@cues <- cbind(CNOlist2@cues, 0)
            CNOlist2@cues <- rbind(CNOlist2@cues,
                                   matrix(0, nrow =nrow(inhibit),
                                          ncol = ncol(CNOlist2@cues)))
            CNOlist2@cues[(nrow(CNOlist2@cues) -
                           nrow(inhibit) + 1):nrow(CNOlist2@cues),
                          which(colnames(CNOlist2@cues) %in%
                                colnames(inhibit))] <- inhibit
            CNOlist2@cues[(nrow(CNOlist2@cues) -
                           nrow(inhibit) + 1):nrow(CNOlist2@cues),
                          ncol(CNOlist2@cues)] <- 1
            colnames(CNOlist2@cues)[ncol(CNOlist2@cues)] <- s
            CNOlist2@stimuli <-
                CNOlist2@cues[, which(colnames(CNOlist2@cues) %in%
                                      colnames(CNOlist2@stimuli))]
            CNOlist2@inhibitors <-
                CNOlist2@cues[, which(colnames(CNOlist2@cues) %in%
                                      c(colnames(CNOlist2@inhibitors), s))]
            CNOlist2@signals[[1]] <-
                rbind(CNOlist2@signals[[1]],
                      matrix(0, nrow = nrow(inhibit),
                             ncol = ncol(CNOlist2@signals[[1]])))
            CNOlist2@signals[[2]] <-
                rbind(CNOlist2@signals[[2]],
                      matrix(0, nrow = nrow(inhibit),
                             ncol = ncol(CNOlist2@signals[[2]])))
        }
        CNOlist2 <- checkCNOlist(CNOlist2)
        return(CNOlist2)
    }
#' @noRd
adj2dnf <-
    function(A) {

        dnf <- NULL

        for (i in seq_len(ncol(A))) {
            for (j in seq_len(nrow(A))) {
                if (i %in% j) { next() }
                if (A[i, j] == 1) {
                    dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j],
                                        sep = "="))
                }
                if (A[i, j] == -1) {
                    dnf <- c(dnf, paste("!", colnames(A)[i], "=",
                                        rownames(A)[j], sep = ""))
                }
            }
        }

        dnf <- unique(dnf)

        return(dnf)

    }
#' @noRd
adj2graph <-
    function(adj.matrix) {
        V   <- rownames(adj.matrix)
        edL <- vector("list", length=nrow(adj.matrix))
        names(edL) <- V
        for (i in seq_len(nrow(adj.matrix))) {
            edL[[i]] <- list(edges=which(!adj.matrix[i,]==0),
                             weights=adj.matrix[i,!adj.matrix[i,]==0])
        }
        gR <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
        return(gR)
    }
#' @noRd
checkCNOlist <-
    function(CNOlist) {
        if (dim(CNOlist@stimuli)[2] > 1) {
            CNOlist@stimuli <- CNOlist@stimuli[,
                                               order(colnames(CNOlist@stimuli))]
        }
        if (dim(CNOlist@inhibitors)[2] > 1) {
            CNOlist@inhibitors <-
                CNOlist@inhibitors[, order(colnames(CNOlist@inhibitors))]
        }
        CNOlist@cues <- cbind(CNOlist@stimuli, CNOlist@inhibitors)
        if (ncol(CNOlist@signals[[1]]) > 1) {
            for (i in seq_len(length(CNOlist@signals))) {
                CNOlist@signals[[i]] <-
                    CNOlist@signals[[i]][,
                                         order(colnames(CNOlist@signals[[i]]))]
            }
        }
        rownames(CNOlist@cues) <- rownames(CNOlist@stimuli) <-
            rownames(CNOlist@inhibitors) <-
            unlist(lapply(rownames(CNOlist@cues), function(x) {
                y <- unlist(strsplit(x, "_"))
                y <- paste(c(sort(y[which(y %in% colnames(CNOlist@stimuli))]),
                             sort(y[which(!(y %in%
                                            colnames(CNOlist@stimuli)))])),
                           collapse = "_")
                return(y)
            }))
        CNOlist@variances <- list()
        return(CNOlist)
    }
#' @noRd
checkMethod <-
    function(method) {
        methods2 <- method
        methods <- c("llr", "cosine", "euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski", "spearman", "pearson",
                     "kendall", "mLL", "cp", "none")
        method <- methods[grep(paste(paste("^", method, sep = ""),
                                     collapse = "|"), methods)[1]]
        method  <- unique(c(method, methods2))
        if (!(any(c("llr", "cosine", "euclidean", "maximum", "manhattan",
                    "canberra", "binary", "minkowski", "spearman", "pearson",
                    "kendall", "mLL", "cp", "none") %in% method))) {
            stop(paste0("I'm sorry, Dave. I'm afraid I can't do that. ",
                        "You have to pick a valid method: ",
                        "llr, cosine, euclidean, maximum, manhattan, ",
                        "canberra, binary, ",
                        "minkowski, spearman, pearson, kendall"))
            ## , mLL, cp, none"))
        }
        return(method)
    }
#' @noRd
#' @importFrom matrixStats rowMins
checkNEMlist <-
    function(NEMlist, CNOlist, parameters, approach, method) {
        NEMlistTmp <- NEMlist
        if("abs" %in% approach) {
            if (length(table(NEMlist$exprs)) == 2) {
                NEMlist$norm <- NEMlist$exprs
                NEMlist$norm[
                            which(NEMlist$norm
                                  ==
                                  as.numeric(names(
                                      table(NEMlist$exprs))[1]))] <-
                    0
                NEMlist$norm[
                            which(NEMlist$norm
                                  ==
                                  as.numeric(
                                      names(table(NEMlist$exprs))[2]))] <-
                    1
            }
            if (length(NEMlist$norm) == 0 &
                !any(c("pearson", "spearman", "kendall") %in% method)) {
                if ("pam" %in% approach) {
                    print("data is not discretized/normed to (0,1);
performing two cluster pam normalization")
                    NEMlist$norm <- pamNorm(NEMlist$exprs)
                } else {
                    print("data is not discretized/normed to (0,1);
performing simple normalization")
                    NEMlist$norm <- simpleNorm(NEMlist$exprs)
                }
                if ("kmeans" %in% approach) {
                    print("data is not discretized/normed to (0,1);
performing two means normalization")
                    NEMlist$norm <- pamNorm(NEMlist$exprs)
                } else {
                    print("data is not discretized/normed to (0,1);
performing simple normalization")
                    NEMlist$norm <- simpleNorm(NEMlist$exprs)
                }
            }
            if ("spearman" %in% method) {
                colnames.exprs <- colnames(NEMlist$exprs)
                rownames.exprs <- rownames(NEMlist$exprs)
                NEMlist$exprs <- t(apply(NEMlist$exprs, 1, rank))
                colnames(NEMlist$exprs) <- colnames.exprs
                rownames(NEMlist$exprs) <- rownames.exprs
            }
        }
        if ("fc" %in% approach) {
            if (length(NEMlist$fc) == 0)  {
                print("foldchanges missing; automatic calculation")
                NEMlist$fc <- computeFc(CNOlist, NEMlist$exprs)
                egenes <- NEMlist$egenes
                if (parameters$cutOffs[2] != 0) {
                    NEMlist <- computeSm(CNOlist, NEMlist, parameters,
                                         method = method)
                }
                NEMlist$egenes <- egenes
            } else {
                ## put somethign here to sort the contrasts abc vs abcd etc...
            }
            if (length(NEMlist$E0) == 0)  {
                egenes <- NEMlist$egenes
                if (parameters$cutOffs[2] != 0) {
                    NEMlist <- computeSm(CNOlist, NEMlist, parameters,
                                         method = method)
                }
                NEMlist$egenes <- egenes
            }
            if ("spearman" %in% method) {
                colnames.fc <- colnames(NEMlist$fc)
                rownames.fc <- rownames(NEMlist$fc)
                NEMlist$fc <- t(apply(NEMlist$fc, 1, rank))
                colnames(NEMlist$fc) <- colnames.fc
                rownames(NEMlist$fc) <- rownames.fc
            }
        }
        if (!is.null(NEMlist$egenes) & is.null(NEMlist$geneGrid)) {
            sgeneAdd <- matrix(Inf, nrow = nrow(NEMlist$fc),
                               ncol = (ncol(CNOlist@signals[[1]])*2))
            for (i in seq_len(nrow(NEMlist$fc))) {
                egeneName <- rownames(NEMlist$fc)[i]
                sgeneCols <- numeric()
                for (j in seq_len(length(NEMlist$egenes))) {
                    if (egeneName %in% NEMlist$egenes[[j]]) {
                        colTmp <- which(colnames(CNOlist@signals[[1]]) ==
                                        names(NEMlist$egenes)[j])
                        sgeneCols <- c(sgeneCols, colTmp,
                        (colTmp+ncol(CNOlist@signals[[1]])))
                    }
                }
                sgeneAdd[i, sgeneCols] <- 0
            }
            sgeneAddCheck <- rowMins(sgeneAdd)
            sgeneAdd[which(sgeneAddCheck == Inf), ] <- 0
            NEMlist$geneGrid <- sgeneAdd
        }
        if (!is.null(NEMlistTmp$weights)) {
            NEMlist$weights <- NEMlistTmp$weights
        }
        NEMlist$signalStates <- NEMlistTmp$signalStates
        return(NEMlist)
    }
#' @noRd
computeFcII <-
    function (y, stimuli = NULL, inhibitors = NULL, batches, runs, extra = 0) {
        design <- makeDesign(y, stimuli, inhibitors, c(batches, runs))
        CompMat <- numeric()
        CompMatNames <- character()
        ctrlsSum <- apply(design[, -grep(paste(c(runs, batches),
                                               collapse = "|"),
                                         colnames(design))], 1, sum)
        ctrlsSum <- which(ctrlsSum == 0)
        stimuliDesign <- design[, grep(paste(stimuli, collapse = "|"),
                                       colnames(design))]
        inhibitorsDesign <- design[, grep(paste(inhibitors, collapse = "|"),
                                          colnames(design))]
        if (length(stimuli) == 1) {
            stimuliDesign <-
                as.matrix(design[, grep(paste(stimuli, collapse = "|"),
                                        colnames(design))])
            colnames(stimuliDesign) <- stimuli
        }
        if (length(inhibitors) == 1) {
            inhibitorsDesign <-
                as.matrix(design[, grep(paste(inhibitors, collapse = "|"),
                                        colnames(design))])
            colnames(inhibitorsDesign) <- inhibitors
        }
        if (is.null(stimuli) == TRUE) {
            stimuliSum <- numeric(nrow(design))
        } else {
            if (length(stimuli) == 1) {
                stimuliSum <- design[, grep(paste(stimuli, collapse = "|"),
                                            colnames(design))]
            } else {
                stimuliSum <- apply(design[, grep(paste(stimuli,
                                                        collapse = "|"),
                                                  colnames(design))], 1, sum)
            }
        }
        if (is.null(inhibitors) == TRUE) {
            inhibitorsSum <- numeric(nrow(design))
        } else {
            if (length(inhibitors) == 1) {
                inhibitorsSum <- design[, grep(paste(inhibitors,
                                                     collapse = "|"),
                                               colnames(design))]
            } else {
                inhibitorsSum <-
                    apply(design[, grep(paste(inhibitors, collapse = "|"),
                                        colnames(design))], 1, sum)
            }
        }
        cuesSum <-
            apply(design[, grep(paste(c(stimuli, inhibitors), collapse = "|"),
                                colnames(design))], 1, sum)
        maxStim <- max(stimuliSum)
        maxKd <- max(inhibitorsSum)
        for (run in runs) {
            for (batch in batches) {
                targetRows <- intersect(which(design[, run] == 1),
                                        which(design[, batch] == 1))
                if (length(targetRows) == 0) { next() }
                grepCtrl <- intersect(which(cuesSum == 0), targetRows)[1]
                grepStims <- intersect(intersect(which(stimuliSum != 0),
                                                 which(inhibitorsSum == 0)),
                                       targetRows)
                grepKds <- intersect(intersect(which(stimuliSum == 0),
                                               which(inhibitorsSum != 0)),
                                     targetRows)
                grepStimsKds <- intersect(intersect(which(stimuliSum != 0),
                                                    which(inhibitorsSum != 0)),
                                          targetRows)
                grepStimsStims <- intersect(intersect(which(stimuliSum >= 2),
                                                      which(inhibitorsSum ==
                                                            0)),
                                            targetRows)
                                        # get ctrl_vs_stim:
                for (i in grepStims) {
                    if (is.na(grepCtrl)) { next() }
                    stimNames <-
                        paste(c("Ctrl", "vs",
                                sort(names(which(stimuliDesign[i, ] >= 1))),
                                batch, run), collapse = "_")
                    if (stimNames %in% CompMatNames) {
                        next()
                    } else {
                        CompMatNames <- c(CompMatNames, stimNames)
                        CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
                    }
                }
                                        # get ctrl_vs_kd:
                for (i in grepKds) {
                    if (is.na(grepCtrl)) { next() }
                    kdNames <-
                        paste(c("Ctrl", "vs",
                                sort(names(which(inhibitorsDesign[i, ] >= 1))),
                                batch, run), collapse = "_")
                    if (kdNames %in% CompMatNames) {
                        next()
                    } else {
                        CompMatNames <- c(CompMatNames, kdNames)
                        CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
                    }
                }
                                        # get kd_vs_kd_stim:
                for (i in grepKds) {
                    kdNames <-
                        paste(sort(names(which(inhibitorsDesign[i, ] >= 1))),
                              collapse = "_")
                    for (j in grepStimsKds) {
                        if (paste(sort(names(which(inhibitorsDesign[j, ] >=
                                                   1))),
                                  collapse = "_") %in% kdNames) {
                            stimNames <-
                                paste(c(kdNames, "vs", kdNames,
                                        sort(names(
                                            which(stimuliDesign[j, ] >= 1))),
                                        batch, run), collapse = "_")
                            if (stimNames %in% CompMatNames) {
                                next()
                            } else {
                                CompMatNames <- c(CompMatNames, stimNames)
                                CompMat <- cbind(CompMat, (y[, j] - y[, i]))
                            }
                        } else {
                            next()
                        }
                    }
                }
                                        # get stim_vs_stim_kd:
                for (i in grepStims) {
                    if (is.na(grepCtrl)) { next() }
                    stimNames <-
                        paste(sort(names(which(stimuliDesign[i, ] >= 1))),
                              collapse = "_")
                    for (j in grepStimsKds) {
                        if (paste(sort(names(which(stimuliDesign[j, ] >= 1))),
                                  collapse = "_") %in% stimNames) {
                            kdNames <-
                                paste(c(stimNames, "vs", stimNames,
                                        sort(names(
                                            which(inhibitorsDesign[j, ] >= 1))),
                                        batch, run), collapse = "_")
                            if (kdNames %in% CompMatNames) {
                                next()
                            } else {
                                CompMatNames <- c(CompMatNames, kdNames)
                                CompMat <- cbind(CompMat, (y[, j] - y[, i]))
                            }
                        } else {
                            next()
                        }
                    }
                }
                if (extra == 1) {
                    ## get stim_vs_stim_stim:
                    for (i in grepStims) {
                        if (is.na(grepCtrl)) { next() }
                        stimNames <-
                            paste(sort(names(which(stimuliDesign[i, ] >= 1))),
                                  collapse = "_")
                        for (j in grepStimsStims) {
                            if (length(unlist(strsplit(stimNames, "_")) %in%
                                       names(which(stimuliDesign[j, ] >= 1)))
                                == length(unlist(strsplit(stimNames, "_"))) &
                                   length(unlist(strsplit(stimNames, "_")) %in%
                                          names(which(stimuliDesign[j, ] >= 1)))
                                < length(names(which(stimuliDesign[j, ] >= 1))))
                            {
                                stim2Names <-
                                    paste(c(stimNames, "vs",
                                            sort(names(which(stimuliDesign[j, ]
                                                             >= 1))),
                                            batch, run), collapse = "_")
                                if (stim2Names %in% CompMatNames) {
                                    next()
                                } else {
                                    CompMatNames <- c(CompMatNames, stim2Names)
                                    CompMat <- cbind(CompMat, (y[, j] - y[, i]))
                                }
                            } else {
                                next()
                            }
                        }
                    }
                    ## get ctrl_vs_stim_kd:
                    for (i in grepStimsKds) {
                        if (is.na(grepCtrl)) { next() }
                        stimNames <-
                            paste(sort(names(which(stimuliDesign[i, ] >= 1))),
                                  collapse = "_")
                        stimkdNames <-
                            paste(c("Ctrl_vs", stimNames,
                                    sort(names(
                                        which(inhibitorsDesign[i, ] >= 1))),
                                    batch, run), collapse = "_")
                        if (stimkdNames %in% CompMatNames) {
                            next()
                        } else {
                            CompMatNames <- c(CompMatNames, stimkdNames)
                            CompMat <- cbind(CompMat, (y[, i] - y[, grepCtrl]))
                        }
                    }
                }
                                        # get s - sk - (k - ctrl): not trivial
                ## for (i in stimuli) {
                ##   for (j in inhibitors) {
                ##     name <-
                ##       CompMatNames <- c(CompMatNames, )
                ##   }
                ## }
            }
        }
        colnames(CompMat) <- CompMatNames
        CompMatCont <- CompMat[, sort(colnames(CompMat))]
        CompMat <- CompMatCont
        return(CompMat)
    }
#' @noRd
#' @import CellNOptR
computeScoreNemT1 <-
    function(CNOlist,
             model,
             bString,
             sizeFac = 10^-10,
             NAFac = 1,
             parameters = list(cutOffs = c(0,1,0), scoring = c(0.25,0.5,2)),
             approach = "fc",
             NEMlist = NULL,
             tellme = 0,
             sim = 0,
             relFit = FALSE,
             method = "s",
             verbose = TRUE,
             opt = "min"
             ) {
        cutModel2 <- function (model, bString) {
            if (sum(bString == 1) > 0) {
                bs = as.logical(bString)
                newmodel <- list()
                if (is.null(dim(model$interMat))) {
                    newmodel$interMat <- model$interMat[bs]
                    newmodel$notMat <- model$notMat[bs]
                    newmodel$reacID <- model$reacID[bs]
                    newmodel$namesSpecies <- model$namesSpecies
                } else {
                    newmodel$interMat <- model$interMat[, bs]
                    newmodel$notMat <- model$notMat[, bs]
                    newmodel$reacID <- model$reacID[bs]
                    newmodel$namesSpecies <- model$namesSpecies
                }
            } else {
                newmodel <- model
            }
            return(newmodel)
        }

        method <- checkMethod(method)
        ## if (is.null(simList) == TRUE) {
        ##   simList = prep4sim(model)
        ## }
        ## if (is.null(indexList) == TRUE) {
        ##   indexList = indexFinder(CNOlist, model)
        ## }
        if (sim == 0) {
            timeIndex = 2
            modelCut = cutModel2(model, bString)
            if (all(bString == 0)) {
                simResults <-
                    simulateStatesRecursive(CNOlist = CNOlist,
                                            model = model, bString = bString,
                                            NEMlist)
            } else {
                simResults <-
                    simulateStatesRecursive(CNOlist = CNOlist, model = modelCut,
                                            bString =
                                                (numeric(
                                                length(modelCut$reacID)) + 1),
                                            NEMlist)
            }
        }
        nInTot <- length(unlist(strsplit(model$reacID, "\\+")))
        nTmp <- unlist(strsplit(model$reacID[which(bString == 1)], "\\+"))
        nInputsNeg <- length(grep("!", nTmp))
        nInputs <- length(nTmp) - nInputsNeg
        ## nInTot <- 1 # use this for no normalization
        sizePen <- sum((c(nInputs, nInputsNeg)/nInTot)*sizeFac)
                                        #/nrow(NEMlist$fc)
        suppressWarnings(Score <- getNemFit(simResults = simResults,
                                            CNOlist = CNOlist, model = modelCut,
                                            timePoint = timeIndex,
                                            NAFac = NAFac, sizePen = sizePen,
                                            simResultsT0 = NA,
                                            approach = approach,
                                            NEMlist = NEMlist, tellme = tellme,
                                            parameters = parameters, sim = sim,
                                            relFit = relFit, method = method,
                                            verbose = verbose, opt = opt))
        if (!is(CNOlist, "CNOlist")) {
            CNOlist = CNOlist(CNOlist)
        }
        return(Score)
    }
#' @noRd
computeSm <-
    function (CNOlist,
              NEMlist,
              parameters,
              method="standard"
              ) {
        CompMatCont <- NEMlist$fc
        fitScore <- -1 # match of high degree
        errorScore <- parameters$scoring[3] # mismatch of high degree
        fitMult <- parameters$scoring[2] # multiplicator for match of low degree
        errorMult <- parameters$scoring[2] # multi. mismatch of low degree
        zeroMult <- parameters$scoring[1]

        CompMat <- CompMatCont
        Epos <- CompMat
        Eneg <- CompMat
        E0 <- CompMat
        EposI <- CompMat*(-1)
        EnegI <- CompMat*(-1)

        beta <- parameters$cutOffs[2]
        alpha <- parameters$cutOffs[1]

        wrongPos <- which(Epos < -beta)
        rightPosII <- which(Epos > alpha & Epos < beta)
        rightPosI <- which(Epos >= beta)
        zerosPosI <- which(Epos <= alpha & Epos >= -alpha)
        zerosPosII <- which(Epos < -alpha & Epos > -beta)

        wrongNeg <- which(Eneg > beta)
        rightNegII <- zerosPosII # which(Eneg < -alpha & Eneg > -beta)
        rightNegI <- which(Eneg <= -beta)
        zerosNegI <- zerosPosI # which(Eneg >= -alpha & Eneg <= alpha)
        zerosNegII <- rightPosII # which(Eneg > alpha & Eneg < beta)

        right0I <- which(abs(E0) <= alpha)
        right0II <- which(abs(E0) > alpha & abs(E0) < beta)
        wrong0 <- which(abs(E0) >= beta)

        if ("mLL" %in% method | "cp" %in% method) {

            E0 <- E0*0
            E0[right0I] <- 1

            Epos <- Epos*0
            Epos[rightPosI] <- 1

            Eneg <- Eneg*0
            Eneg[rightNegI] <- 1

            EposI <- EposI*0
            EposI[rightPosII] <- 1

            EnegI <- EnegI*0
            EnegI[rightNegII] <- 1

        } else {

            E0[wrong0] <- errorScore*zeroMult
            E0[right0II] <- fitScore*zeroMult*fitMult
            E0[right0I] <- fitScore*zeroMult

            Epos[zerosPosI] <- errorScore*errorMult*zeroMult
            Epos[zerosPosII] <- errorScore*errorMult
            Epos[wrongPos] <- errorScore
            Epos[rightPosI] <- fitScore
            Epos[rightPosII] <- fitScore*fitMult

            Eneg[zerosNegI] <- -errorScore*errorMult*zeroMult
            Eneg[zerosNegII] <- -errorScore*errorMult
            Eneg[wrongNeg] <- -errorScore
            Eneg[rightNegI] <- -fitScore
            Eneg[rightNegII] <- -fitScore*fitMult

            EposI[zerosNegI] <- errorScore*errorMult*zeroMult
            EposI[zerosNegII] <- errorScore*errorMult
            EposI[wrongNeg] <- errorScore
            EposI[rightNegI] <- fitScore
            EposI[rightNegII] <- fitScore*fitMult

            EnegI[zerosPosI] <- -errorScore*errorMult*zeroMult
            EnegI[zerosPosII] <- -errorScore*errorMult
            EnegI[wrongPos] <- -errorScore
            EnegI[rightPosI] <- -fitScore
            EnegI[rightPosII] <- -fitScore*fitMult

        }

        return(list(exprs = NEMlist$exprs, fc = CompMatCont, E0 = E0,
                    Epos = Epos, Eneg = Eneg, EposI = EposI, EnegI = EnegI))
    }
#' @noRd
deleteEdge <-
    function(edges, CNOlist, model, n = 100) {
        sifMatrix <- numeric()
        graph <- model$reacID[-grep("\\+", model$reacID)]
        for (i in edges) {
            if (i %in% graph) {
                graph <- graph[-which(graph %in% i)]
            }
        }
        for (i in graph) {
            tmp2 <- unlist(strsplit(i, "="))
            if (gsub("!", "", tmp2[1]) %in% tmp2[1]) {
                sifMatrix <- rbind(sifMatrix, c(tmp2[1], 1, tmp2[2]))
            } else {
                sifMatrix <- rbind(sifMatrix, c(gsub("!", "", tmp2[1]), -1,
                                                tmp2[2]))
            }
        }
        write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
        PKN2 <- readSIF("temp.sif")
        unlink("temp.sif")
        checkSignals(CNOlist,PKN2)
        indices<-indexFinder(CNOlist,PKN2,verbose=FALSE)
        NCNOindices<-findNONC(PKN2,indices,verbose=FALSE)
        NCNOcut<-cutNONC(PKN2,NCNOindices)
        indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
        NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
        indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
        NCNOcutCompExp<-expandNEM(NCNOcutComp, maxInputsPerGate = n)
        resECNOlist<-residualError(CNOlist)
        Fields4Sim<-prep4sim(NCNOcutCompExp)
        model2 <- NCNOcutCompExp
        return(model2)
    }
#' @noRd
deleteSignal <-
    function(s, CNOlist) {
        CNOlist2 <- CNOlist
        for (i in s) {
            if (i %in% colnames(CNOlist2@signals[[1]])) {
                CNOlist2@signals[[1]] <-
                    CNOlist2@signals[[1]][
                               , -which(colnames(CNOlist2@signals[[1]]) %in% i)]
                CNOlist2@signals[[2]] <-
                    CNOlist2@signals[[2]][
                               , -which(colnames(CNOlist2@signals[[2]]) %in% i)]
            }
            if (i %in% colnames(CNOlist2@stimuli)) {
                CNOlist2@cues <-
                    CNOlist2@cues[, -which(colnames(CNOlist2@cues) %in% i)]
                if (!is.null(
                         dim(CNOlist2@stimuli[
                                        , -which(colnames(CNOlist2@stimuli)
                                                 %in% i)]))) {
                    CNOlist2@stimuli <-
                        CNOlist2@stimuli[
                                   , -which(colnames(CNOlist2@stimuli) %in% i)]
                } else {
                    CNOlist2@stimuli <-
                        as.matrix(CNOlist2@stimuli[
                                             , -which(colnames(CNOlist2@stimuli)
                                                      %in% i)])
                }
            }
            if (i %in% colnames(CNOlist2@inhibitors)) {
                CNOlist2@cues <-
                    CNOlist2@cues[, -which(colnames(CNOlist2@cues) %in% i)]
                if (!is.null(
                         dim(CNOlist2@inhibitors[
                                        , -which(colnames(CNOlist2@inhibitors)
                                                 %in% i)]))) {
                    CNOlist2@inhibitors <-
                        CNOlist2@inhibitors[
                                   , -which(colnames(CNOlist2@inhibitors)
                                            %in% i)]
                } else {
                    CNOlist2@inhibitors <-
                        as.matrix(
                            CNOlist2@inhibitors[
                                       , -which(colnames(CNOlist2@inhibitors)
                                                %in% i)])
                }
            }
        }
        return(CNOlist2)
    }
#' @noRd
disc <-
    function(x, beta = 0.5) { # simple discretization
        x.disc <- x
        x.disc[which(abs(x.disc) <= beta)] <- 0
        x.disc[which(x.disc > beta)] <- 1
        x.disc[which(x.disc < -beta)] <- -1
        return(x.disc)
    }
#' @noRd
dnf2adj <-
    function(dnf, closed = FALSE) {
        if (length(dnf) == 0) {
            return(NULL)
        } else {
            nodes <- character()
            for (i in dnf) {
                tmp <- unlist(strsplit(i, "="))
                nodes <- c(nodes, tmp[2], unlist(strsplit(tmp[1], "\\+")))
            }
            nodes <- unique(gsub("!", "", nodes))
            adjmat <- matrix(0, length(nodes), length(nodes))
            colnames(adjmat) <- nodes
            rownames(adjmat) <- nodes
            for (i in dnf) {
                tmp <- unlist(strsplit(i, "="))
                child <- tmp[2]
                parents <- unlist(strsplit(tmp[1], "\\+"))
                for (j in parents) {
                    if (gsub("!", "", j) %in% j) {
                        adjmat[which(rownames(adjmat) %in% j),
                               which(colnames(adjmat) %in% child)] <- 1
                    } else {
                        adjmat[which(rownames(adjmat) %in% gsub("!", "", j)),
                               which(colnames(adjmat) %in% child)] <- -1
                    }
                }
            }
            diag(adjmat) <- 1
            if (closed) {
                stop <- FALSE
                cons <- c(TRUE, rep(FALSE, (length(adjmat) - 1)))
                while(!stop) {
                    adjmat <- adjmat%*%adjmat
                    if (all(cons == (adjmat != 0))) {
                        stop <- TRUE
                    } else {
                        cons <- (adjmat != 0)
                    }
                }
                adjmat[adjmat > 1] <- 1
                adjmat[adjmat < -1] <- -1
            }
            return(adjmat)
        }
    }
#' @noRd
#' @importFrom Rgraphviz lines
drawLocal <-
    function(l, type = "l", ...) {
        local.edges <- c(1, l$edges[[1]])
        local.edges1 <- local.edges
        if (length(unique(local.edges)) < 2) {
            local.edges <- rep(mean(l$scores[[1]]), length(local.edges1))
            min.local.edges <- min(local.edges)
            convert.edges <- seq(min(l$scores[[1]]), max(l$scores[[1]]),
                                 length.out = max(local.edges -
                                                  min.local.edges))
            convert.legend <- sort(unique(local.edges))
        } else {
            min.local.edges <- min(local.edges)
            convert.edges <- seq(min(l$scores[[1]]), max(l$scores[[1]]),
                                 length.out = max(local.edges -
                                                  min.local.edges)+1)
            convert.legend <- sort(unique(local.edges))
            for (i in seq_len((max(local.edges - min.local.edges)+1))) {
                local.edges[which(local.edges == (i+min.local.edges)-1)] <-
                    convert.edges[which(convert.legend ==
                                        (i+min.local.edges)-1)]
            }
        }
        ylim <- c(min(l$scores[[1]]), max(l$scores[[1]]))
        par(mar=c(4, 4, 4, 4) + 0.1)
        plot(l$scores[[1]], type = type, main = "score improvement (black)
and number of changed edges (red)", ylab = "score", xlab = "move", ...)
        lines(local.edges, col = 2, type = "b")
        axis(4, at = sort(unique(local.edges)),
             labels = sort(unique(local.edges1)), col = "red",
             col.ticks = "red", col.axis = "red")
        mtext("edges changed", side=4, line=3, cex.lab=1,las=0, col="red")
    }
#' @noRd
#' @importFrom Rgraphviz lines
drawScores <-
    function(CNOresult) {
        gens <- CNOresult$results[, 1]
        stallGens <- as.numeric(unlist(strsplit(gens, " / ")))
        stallGens <- stallGens[(seq_len((length(stallGens)/2)))*2 - 1]
        getInfo <- CNOresult$results[which(stallGens == 0), ]
        bestScore <- unlist(strsplit(getInfo[, 2], " "))
        bestScore <- as.numeric(gsub("\\(", "", gsub("\\)", "", bestScore)))
        bestScore <- bestScore[(seq_len((length(bestScore)/2)))*2 - 1]
        avgScore <- unlist(strsplit(getInfo[, 2], " "))
        avgScore <- as.numeric(
            gsub("\\)", "", gsub("\\(", "",
                                 avgScore[(seq_len((length(avgScore)/2)))*2])))
        par(mfrow=c(3,1))
        plot(seq_len(length(bestScore)), bestScore, col = "red", type = "l",
             main = paste("Score Improvement", sep = ""), xlab = "Generation",
             ylab = "Score",
             ylim = c(min(bestScore), max(c(bestScore, avgScore))), xaxt='n')
        axis(1, at = seq_len(length(bestScore)), labels = which(stallGens == 0))
        lines(seq_len(length(bestScore)), avgScore, col = "blue", type = "l")
        legend(length(bestScore), max(c(bestScore, avgScore)),
               legend = c("Best Score", "Average Score"),
               fill = c("red", "blue"), xjust = 1)
        if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 200) {
            overallTime <-
                paste(
                    round(
                        sum(as.numeric(
                            CNOresult$results[, "Iter_time"]))/60, 3),
                    " minutes", sep = "")
        }
        if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 18000) {
            overallTime <-
                paste(
                    round(
                        sum(as.numeric(
                            CNOresult$results[, "Iter_time"]))/60/60, 4),
                    " hours", sep = "")
        }
        if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 864000) {
            overallTime <-
                paste(
                    round(
                        sum(as.numeric(
                            CNOresult$results[, "Iter_time"]))/60/60/24, 5),
                    " days", sep = "")
        }
        if (sum(as.numeric(CNOresult$results[, "Iter_time"])) <= 200 ) {
            overallTime <-
                paste(
                    round(
                        sum(as.numeric(CNOresult$results[, "Iter_time"])), 2),
                    " seconds", sep = "")
        }
        plot(as.numeric(CNOresult$results[, "Iter_time"]),
             xlab = "Generation", ylab = "Iteration Time in seconds",
             main = paste("Iteration Time over all Generations (time overall: ",
                          overallTime, ")", sep = ""), type = "l")
        bStringDist <- numeric()
        for (i in seq_len((nrow(getInfo) - 1))) {
            bStringDist <-
                c(bStringDist,
                  dist(rbind(unlist(strsplit(getInfo[i, 3], ", ")),
                             unlist(strsplit(getInfo[i+1, 3], ", ")))))
        }
        bStringDist <- c(0, bStringDist)
        plot(seq_len(length(bStringDist)), bStringDist, type = "l",
             main = paste("Best Strings Distance (average distance: ",
                          round(mean(bStringDist), 2), ")", sep = ""),
             xlab = "Generation", ylab = "Euclidean Distance", xaxt='n')
        axis(1, at = seq_len(length(bestScore)), labels = which(stallGens == 0))
        diffs <- numeric()
        for (i in seq_len(length(unlist(strsplit(getInfo[1, 3], ", "))))) {
            diffs <- c(diffs, dist(rbind(rep(0, i), rep(1, i))))
        }
        axis(4, at = diffs, labels = seq_len(length(diffs)))
    }
#' @noRd
expandNEM <-
    function(model, ignoreList=NA,maxInputsPerGate=2){

        Model = model
                                        # check that Model is a Model list
        if(!is.list(Model)) stop("This function expects as input a Model
as output by readSIF")
        if(length(Model) == 4) {
            if(all(names(Model) != c("reacID",
                                     "namesSpecies","interMat","notMat"))){
                stop("This function expects as input a Model
as output by readSIF")
            }
        }
        if(length(Model) == 5) {
            if(all(names(Model) != c("reacID","namesSpecies","interMat",
                                     "notMat","speciesCompressed"))) {
                stop("This function expects as input a Model
as output by readSIF")
            }
        }

        SplitANDs <- list(initialReac=c("split1","split2"))
        splitR <- 1

        ## split all the ANDs
        ## remove any ANDs 2/3 and save >3 to add later
        ## +3 won't get added here again but are part of prior
        ## knowledge if contained within PKN
        andToAdd = c()
        remove.and = c()
        reacs2Ignore = c()
        initialReacN <- length(Model$reacID)

        ## TODO: move into readSIF ?
        if (initialReacN == 1){
            Model$interMat <- as.matrix(Model$interMat)
        }

        ## which reactions have ignoreList as output?
        if(!is.na(ignoreList[1])) {
            for(s in seq_len(initialReacN)) {
                if(any(Model$interMat[ignoreList,s] == 1)) {
                    reacs2Ignore = c(reacs2Ignore, s)
                }
            }
        }

        for(r in seq_len(initialReacN)) {
            inNodes <- which(Model$interMat[,r] == -1)
            if(length(inNodes) > 1) {
                if(length(inNodes) > 3) {
                    andToAdd = c(andToAdd, r)
                }
                remove.and = c(remove.and, r)

                if(!any(reacs2Ignore == r)) {
                    outNode <- which(Model$interMat[,r] == 1)
                    newReacs <- matrix(data=0,nrow=dim(Model$interMat)[1],
                                       ncol=length(inNodes))
                    newReacsNots <- matrix(data=0,nrow=dim(Model$interMat)[1],
                                           ncol=length(inNodes))
                    newReacs[outNode,] <- 1
                    newReacsIDs <- rep("a",length(inNodes))
                    for(i in seq_len(length(inNodes))) {
                        newReacs[inNodes[i],i] <- -1
                        newReacsNots[inNodes[i],i] <- Model$notMat[inNodes[i],r]
                        newReacsIDs[i] <-
                            paste(Model$namesSpecies[inNodes[i]],"=",
                                  Model$namesSpecies[outNode],sep="")
                        if (Model$notMat[inNodes[i],r] == 1) {
                            newReacsIDs[i] <- paste("!",newReacsIDs[i],sep="")
                        }
                    }
                    colnames(newReacs) <- newReacsIDs
                    colnames(newReacsNots) <- newReacsIDs
                    SplitANDs[[splitR]] <- newReacsIDs
                    names(SplitANDs)[splitR] <- Model$reacID[r]
                    splitR <- splitR+1
                    Model$notMat <- cbind(Model$notMat,newReacsNots)
                    Model$interMat <- cbind(Model$interMat,newReacs)
                    Model$reacID <- c(Model$reacID,newReacsIDs)
                }
            }
        }

        if(length(andToAdd)) {
            toAdd = list()
            toAdd$notMat <- Model$notMat[,andToAdd]
            toAdd$interMat <- Model$interMat[,andToAdd]
            toAdd$reacID <- Model$reacID[andToAdd]
        } else {
            toAdd <- NA
        }

        if(length(remove.and)) {
            Model$notMat <- Model$notMat[,-remove.and]
            Model$interMat <- Model$interMat[,-remove.and]
            Model$reacID <- Model$reacID[-remove.and]
        }

        newANDs <- list(finalReac=c("or1","or2"))
        ANDsadded <- 1
        total.list = seq_len(length(Model$namesSpecies))


        ## functions to get lhs and rhs of reactions
        getlhs <- function(x) {
            spec1 = strsplit(x, "=")[[1]][1]
        }
        getrhs <- function(x) {
            spec2 = strsplit(x, "=")[[1]][2]
        }

        ## scan all species and build and gates if required
        for(sp in total.list) {
            inReacsIndex <- which(Model$interMat[sp,] == 1)
            if(length(inReacsIndex) > 1) {
                inReacs <- Model$interMat[,inReacsIndex]
                findInput <- function(x) {
                    inNode<-which(x == -1)
                }
                inSp <- apply(inReacs,2,findInput)

                ## let
                ## find the input species first and store in a vector
                inSpecies = apply(as.matrix(colnames(inReacs)), 1, getlhs)
                outname = Model$namesSpecies[sp]

                ## just for sanity check, all outputs must be the same
                outnames = apply(as.matrix(colnames(inReacs)), 1, getrhs)
                if (length(unique(outnames))!=1 | outname!=outnames[1]){
                    stop("error in expandGates.
should not happen here. please report")
                }

                                        # an alias
                myrownames = rownames(Model$interMat)

                                        # first the 2 inputs cases

                combinations = combn(seq(1,length(inSpecies)), 2)

                for (this in seq(1, dim(combinations)[2])){
                    i = combinations[1,this]
                    j = combinations[2,this]
                    ## names[i] and names[j] contains possibly the !
                    ## sign,let us get the real species names
                    realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",
                                       substr(inSpecies[i],2,10000),
                                       inSpecies[i])
                    realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",
                                       substr(inSpecies[j],2,10000),
                                       inSpecies[j])

                    realnames = c(realname1,realname2)
                    if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){
                        ## exclude reaction if any name are indentical
                        next()
                    }

                    ## create the new reaction Id to be used as a column name
                    newcolname = paste(paste(inSpecies[i],
                                             inSpecies[j],sep="+"),
                                       outname, sep="=")
                    if (newcolname %in% colnames(Model$interMat)){
                        next() # skip if exist already
                    }
                    Model$reacID <- c(Model$reacID,newcolname)

                    ## fill the interMat (-1 if in lhs, 1 if in rhs)
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    values[which(myrownames == realname1)]<- -1
                    values[which(myrownames == realname2)]<- -1
                    values[which(myrownames == outname)]<- 1
                    Model$interMat= cbind(Model$interMat, values)

                                        # now, the notMat, 0 by default
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    if (substr(inSpecies[i],1,1) == "!"){ #look first specy
                        values[which(myrownames == realname1)]<- 1
                    }
                    if (substr(inSpecies[j],1,1) == "!"){# and second one
                        values[which(myrownames == realname2)]<- 1
                    }
                    Model$notMat= cbind(Model$notMat, values)

                                        # finally, fill the newAnd list
                    newreac1 = paste(inSpecies[i], outname, sep="=")
                    newreac2 = paste(inSpecies[j], outname, sep="=")
                    newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2)
                    names(newANDs)[[length(newANDs)]] <- newcolname
                }

                ## Same code as above but to create the 3 inputs combinations
                if (length(inSpecies)>=3 & maxInputsPerGate>=3){
                    combinations = combn(seq(1,length(inSpecies)), 3)
                    indices = seq(1,dim(combinations)[2])
                }
                else{
                    indices = seq(length=0)
                }

                for (this in indices){
                    i = combinations[1,this]
                    j = combinations[2,this]
                    k = combinations[3,this]
                    realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",
                                       substr(inSpecies[i],2,10000),
                                       inSpecies[i])
                    realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",
                                       substr(inSpecies[j],2,10000),
                                       inSpecies[j])
                    realname3 = ifelse(substr(inSpecies[k], 1,1) =="!",
                                       substr(inSpecies[k],2,10000),
                                       inSpecies[k])

                    realnames = c(realname1,realname2, realname3)
                    if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){
                        ## exclude reaction if any name are indentical
                        next()
                    }
                    newcolname <- paste(paste(inSpecies[i],
                                              inSpecies[j],inSpecies[k],
                                              sep="+"), outname, sep="=")
                    if (newcolname %in% colnames(Model$interMat)){
                        next() # skip if exist already
                    }
                    Model$reacID <- c(Model$reacID,newcolname)

                                        # intermat first
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    for (name in inSpecies){
                        realname = ifelse(substr(name, 1,1)=="!",
                                          substr(name,2,10000), name)
                        values[which(myrownames == realname)]<- -1
                    }
                    values[which(myrownames == outname)]<- 1
                    Model$interMat= cbind(Model$interMat, values)

                                        # now, the notMat
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    for (name in inSpecies){
                        if (substr(name,1,1) == "!"){
                            realname = ifelse(substr(name, 1,1) =="!",
                                              substr(name,2,10000), name)
                            values[which(myrownames == realname)]<- 1
                        }
                    }
                    Model$notMat= cbind(Model$notMat, values)

                                        # finally the newAnd
                    newreac1 = paste(inSpecies[i], outname, sep="=")
                    newreac2 = paste(inSpecies[j], outname, sep="=")
                    newreac3 = paste(inSpecies[k], outname, sep="=")
                    newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2,
                                                      newreac3)
                    names(newANDs)[[length(newANDs)]] <- newcolname
                }

                ## Same code as above but to create the 4 inputs combinations
                if (length(inSpecies)>=4 & maxInputsPerGate>=4){
                    combinations = combn(seq(1,length(inSpecies)), 4)
                    indices = seq(1,dim(combinations)[2])
                }
                else{
                    indices = seq(length=0)
                }

                for (this in indices){
                    i = combinations[1,this]
                    j = combinations[2,this]
                    k = combinations[3,this]
                    l = combinations[4,this]
                    realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",
                                       substr(inSpecies[i],2,10000),
                                       inSpecies[i])
                    realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",
                                       substr(inSpecies[j],2,10000),
                                       inSpecies[j])
                    realname3 = ifelse(substr(inSpecies[k], 1,1) =="!",
                                       substr(inSpecies[k],2,10000),
                                       inSpecies[k])
                    realname4 = ifelse(substr(inSpecies[l], 1,1) =="!",
                                       substr(inSpecies[l],2,10000),
                                       inSpecies[l])
                    realnames = c(realname1,realname2, realname3, realname4)
                    if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){
                        ## exclude reaction if any name are indentical
                        next()
                    }
                    newcolname <- paste(paste(inSpecies[i], inSpecies[j],
                                              inSpecies[k],inSpecies[l],
                                              sep="+"), outname, sep="=")
                    if (newcolname %in% colnames(Model$interMat)){
                        next() # skip if exist already
                    }
                    Model$reacID <- c(Model$reacID,newcolname)

                                        # intermat first
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    for (name in inSpecies){
                        realname = ifelse(substr(name, 1,1)=="!",
                                          substr(name,2,10000), name)
                        values[which(myrownames == realname)]<- -1
                    }
                    values[which(myrownames == outname)]<- 1
                    Model$interMat= cbind(Model$interMat, values)

                                        # now, the notMat
                    values = as.matrix(rep(0, length(Model$namesSpecies)))
                    colnames(values)<-newcolname
                    for (name in inSpecies){
                        if (substr(name,1,1) == "!"){
                            realname = ifelse(substr(name, 1,1) =="!",
                                              substr(name,2,10000), name)
                            values[which(myrownames == realname)]<- 1
                        }
                    }
                    Model$notMat= cbind(Model$notMat, values)

                                        # finally the newAnd
                    newreac1 = paste(inSpecies[i], outname, sep="=")
                    newreac2 = paste(inSpecies[j], outname, sep="=")
                    newreac3 = paste(inSpecies[k], outname, sep="=")
                    newreac4 = paste(inSpecies[l], outname, sep="=")
                    newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2,
                                                      newreac3, newreac4)
                    names(newANDs)[[length(newANDs)]] <- newcolname

                } # end if length(inSp) == 2
                if (maxInputsPerGate >= 5) {
                    for (mip in 5:maxInputsPerGate) {
                        if (length(inSpecies) >= mip & maxInputsPerGate >=
                            mip) {
                            combinations = combn(seq(1, length(inSpecies)),
                                                 mip)
                            indices = seq(1, dim(combinations)[2])
                        }
                        else {
                            indices = seq(length = 0)
                        }
                        for (this in indices) {
                            combs <- c()
                            realnames <- c()
                            for (i in seq_len(mip)) {
                                combs[i] = combinations[i, this]
                                realnames[i] =
                                    ifelse(substr(inSpecies[combs[i]], 1, 1) ==
                                           "!", substr(inSpecies[i],
                                                       2, 10000),
                                           inSpecies[combs[i]])
                            }
                            if (any(combn(realnames, 2)[1, ] ==
                                    combn(realnames,
                                          2)[2, ])) {
                                (next)()
                            }
                            newcolname <- paste(paste(inSpecies[combs],
                                                      collapse = "+"), outname,
                                                sep = "=")
                            if (newcolname %in% colnames(Model$interMat)) {
                                (next)()
                            }
                            Model$reacID <- c(Model$reacID, newcolname)
                            values = as.matrix(rep(0,
                                                   length(Model$namesSpecies)))
                            colnames(values) <- newcolname
                            for (name in inSpecies) {
                                realname = ifelse(substr(name, 1, 1) == "!",
                                                  substr(name, 2, 10000), name)
                                values[which(myrownames == realname)] <- -1
                            }
                            values[which(myrownames == outname)] <- 1
                            Model$interMat = cbind(Model$interMat, values)
                            values = as.matrix(rep(0,
                                                   length(Model$namesSpecies)))
                            colnames(values) <- newcolname
                            for (name in inSpecies) {
                                if (substr(name, 1, 1) == "!") {
                                    realname = ifelse(substr(name, 1, 1) == "!",
                                                      substr(name, 2, 10000),
                                                      name)
                                    values[which(myrownames == realname)] <- 1
                                }
                            }
                            Model$notMat = cbind(Model$notMat, values)
                            newreac1 = paste(inSpecies[i], outname, sep = "=")
                            newreac2 = paste(inSpecies[j], outname, sep = "=")
                            newreac3 = paste(inSpecies[k], outname, sep = "=")
                            newreac4 = paste(inSpecies[l], outname, sep = "=")
                            newreacs <- c()
                            for (i in seq_len(mip)) {
                                newreacs[i] <- paste(inSpecies[combs[i]],
                                                     outname, sep = "=")
                            }
                            newANDs[[length(newANDs) + 1]] <- newreacs
                            names(newANDs)[[length(newANDs)]] <- newcolname
                        }
                    }
                }

            }
        }

        if(!is.na(toAdd)) {
            Model$notMat = cbind(Model$notMat, toAdd$notMat)
            Model$interMat = cbind(Model$interMat, toAdd$interMat)
            Model$reacID = c(Model$reacID, toAdd$reacID)
        }

        modelExp <- Model
        modelExp$SplitANDs <- SplitANDs
        modelExp$newANDs <- newANDs
        return(modelExp)
    }
#' @noRd
expNorm <-
    function(x,  stimuli = NULL, inhibitors = NULL, batches, runs, cutoff) {
        design <- makeDesign(x, stimuli, inhibitors, c(batches, runs))
        normedX <- x*0
        for (run in runs) {
            for (batch in batches) {
                targetRows <- intersect(which(design[, run] == 1),
                                        which(design[, batch] == 1))
                if (length(targetRows) == 0) { next() }
                for (i in seq_len(nrow(normedX))) {
                    if (max(x[i, targetRows]) -
                        min(x[i, targetRows]) >= cutoff) {
                        normedX[i, targetRows] <-
                            simpleNorm(x[i, targetRows])
                    }
                }
            }
        }
        cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors),
                                             collapse = "|"),
                                       colnames(design))], 1, sum)
        grepCtrl <- which(cuesSum == 0)
        for (run in runs) {
            for (batch in batches) {
                targetRows <- intersect(which(design[, run] == 1),
                                        which(design[, batch] == 1))
                if (length(targetRows) == 0) { next() }
                for (i in seq_len(nrow(normedX))) {
                    if (max(x[i, targetRows]) -
                        min(x[i, targetRows]) < cutoff) {
                        normedX[i, targetRows] <- median(normedX[i, grepCtrl])
                    }
                }
            }
        }
        return(normedX)
    }
#' @noRd
#' @importFrom mnem plotDnf
exportVars <-
    function(type = "ga") { # type: ga, loc, ex
        if ("ga" %in% type) {
            return(list("computeScoreNemT1", "simulateStatesRecursiveAdd",
                        "simulateStatesRecursive", "reduceGraph", "getNemFit",
                        "checkCNOlist", "checkNEMlist", "computeFc",
                        "computeSm", "relFit", "sizeFac", "method",
                        "removeCycles", "dnf2adj", "verbose", "getHierarchy",
                        "absorption", "checkMethod"))
        }
        if ("loc" %in% type) {
            return(list("computeScoreNemT1", "simulateStatesRecursiveAdd",
                        "simulateStatesRecursive", "reduceGraph", "getNemFit",
                        "checkCNOlist", "checkNEMlist", "computeFc",
                        "computeSm", "relFit", "sizeFac", "method",
                        "absorptionII", "max.steps", "max.time", "removeCycles",
                        "node", "dnf2adj", "verbose", "absorpII", "draw",
                        "getHierarchy", "absorption", "bitStrings",
                        "checkMethod"))
        }
        if ("ex" %in% type) {
            return(list("computeScoreNemT1", "simulateStatesRecursiveAdd",
                        "simulateStatesRecursive", "reduceGraph", "getNemFit",
                        "checkCNOlist", "checkNEMlist", "computeFc",
                        "computeSm", "relFit", "sizeFac", "method",
                        "removeCycles", "dnf2adj", "verbose", "getHierarchy",
                        "absorption", "checkMethod", "NEMlist"))
        }
    }
#' @noRd
#' @import snowfall CellNOptR
exSearch <-
    function(CNOlist,model,sizeFac=10^-10,NAFac=1,NEMlist,
             parameters=list(cutOffs=c(0,1,0), scoring=c(0.1,0.2,0.9)),
             parallel = NULL, method = "s", relFit = FALSE, verbose = TRUE,
             reduce = TRUE, approach = "fc", ...) {

        cutModel2 <- function (model, bString) {
            if (sum(bString == 1) > 0) {
                bs = as.logical(bString)
                newmodel <- list()
                if (is.null(dim(model$interMat))) {
                    newmodel$interMat <- model$interMat[bs]
                    newmodel$notMat <- model$notMat[bs]
                    newmodel$reacID <- model$reacID[bs]
                    newmodel$namesSpecies <- model$namesSpecies
                } else {
                    newmodel$interMat <- model$interMat[, bs]
                    newmodel$notMat <- model$notMat[, bs]
                    newmodel$reacID <- model$reacID[bs]
                    newmodel$namesSpecies <- model$namesSpecies
                }
            } else {
                newmodel <- model
            }
            return(newmodel)
        }


        bin2dec <- function(x) {
            exp2 <- 2^c((length(x)-1):0)
            y <- exp2%*%x
            return(y)
        }

        dec2bin <- function(x) {
            if (x == 0) {
                y <- 0
            } else {
                tmp <- rev(as.integer(intToBits(x)))
                y <- tmp[which(tmp != 0)[1]:length(tmp)] # fast dec2bin
            }
            return(y)
        }

        dec2binOld <- function(x) {
            if (x == 0) {
                y <- 0
            } else {
                xTmp <- x
                start <- floor(log(xTmp)/log(2))
                y <- numeric(start)
                count <- 0
                for (i in start:0) {
                    count <- count + 1
                    exp <- floor(log(xTmp)/log(2))
                    if (i == exp) {
                        y[count] <- 1
                        xTmp <- xTmp - 2^exp
                    } else {
                        y[count] <- 0
                    }
                }
            }
            return(y)
        }

        if (!is.null(parallel)) {
            if (is.list(parallel)) {
                if (length(parallel[[1]]) != length(parallel[[2]])) {
                    stop("The nodes (second list object in parallel) and
the number of cores used on every node (first list object in parallel)
must be the same.") }
                hosts <- character()
                for (i in seq_len(length(parallel[[1]]))) {
                    hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
                }
                hosts <- as.list(hosts)
                sfInit(parallel=TRUE, socketHosts=hosts)
            } else {
                sfInit(parallel=TRUE, cpus=parallel)
            }
            sfLibrary(bnem)
        }
        spaceExp <- 2^length(model$reacID)
        print(paste("structurally different networks: ", spaceExp, sep = ""))
        ## reduce to equivalent classes by absorption
        bigBang <- function(j) {
            essential <- dec2bin((j-1))
            tmp <- reduceGraph(c(rep(0,
            (length(model$reacID)-length(essential))), essential), model,
            CNOlist)
            return(tmp)
        }
        pop <- matrix(NA, nrow = spaceExp, ncol = length(model$reacID))
        if (!is.null(parallel)) {
            pop <- sfLapply(as.list(seq_len(spaceExp)), bigBang)
        } else {
            pop <- lapply(as.list(seq_len(spaceExp)), bigBang)
            cat("\n")
        }
        pop <- do.call("rbind", pop)
        res <- pop
        res <- apply(res, 1, paste, collapse = "")
        realSpace <- seq_len(spaceExp)
        if (sum(duplicated(res) == TRUE) > 0 & reduce) {
            realSpace <- realSpace[-which(duplicated(res) == TRUE)]
        }
        ## pop <- pop[-which(duplicated(res) == TRUE), ]
        print(paste("reduction to structurally different equivalence classes: ",
                    length(realSpace), sep = ""))
        scores <- numeric(spaceExp)
        getBinScore <- function(j) {
            cat(".")
            scores[j] <- computeScoreNemT1(CNOlist, model = model, pop[j, ],
                                           sizeFac = sizeFac, NAFac = NAFac,
                                           NEMlist = NEMlist,
                                           parameters = parameters,
                                           method = method, relFit = relFit,
                                           approach=approach)
            return(scores[j])
        }
        if (!is.null(parallel)) {
            res <- sfApply(as.matrix(realSpace), 1, getBinScore)
        } else {
            res <- apply(as.matrix(realSpace), 1, getBinScore)
            cat("\n")
        }
        if (sizeFac == 0) {
            getSize <- function(x) {
                dnf <- model$reacID[as.logical(x)]
                nTmp <- unlist(strsplit(dnf, "\\+"))
                return(length(nTmp))
            }
            sizes <- apply(pop[realSpace, ], 1, getSize)
            minscores <- which(res == min(res))
            minscore <- minscores[order(sizes[minscores],
                                        decreasing = FALSE)[1]]
            bString <- (pop[realSpace, ])[minscore, ]
        } else {
            bString <- (pop[realSpace, ])[which(res == min(res))[1], ]
        }
        print(bString)
        dtmRatio <- mean(abs(res - mean(res)))/abs(min(res) - mean(res))
        dtmRatio2 <- 1 - abs(min(res) - mean(res))/abs(min(res) - max(res))
        print("best network found:")
        print(toString(bString))
        print("score:")
        print(min(res))
        print("distance to mean ratio (smaller is better):")
        print(dtmRatio)
        print(dtmRatio2)
                                        #hist(res)
        if (!is.null(parallel)) {
            sfStop()
        }
        plotDnf(model$reacID[as.logical(bString)])
        ## names(res) <- samples
        return(list(bString = bString, score = min(res), bStrings = pop,
                    scores = res, dtmRatio = dtmRatio))
    }
#' @noRd
#' @importFrom mnem plotDnf
#' @importFrom Rgraphviz lines
#' @import snowfall
gaBinaryNemT1 <-
    function (CNOlist,
              model,
              initBstring = NULL, # initBstring = TRUE
              sizeFac = 10^-10,
              NAFac = 1,
              popSize = 100,
              pMutation = 0.5,
              maxTime = Inf,
              maxGens = Inf,
              stallGenMax = 10,
              relTol = 0.01,
              verbose = TRUE,
              priorBitString = NULL,
              selPress = c(1.2,0.0001), # 1.2
              approach = "fc",
              NEMlist,
              fit = "linear",
              targetBstring = "none",
              elitism = NULL,
              inversion = NULL,
              graph = TRUE,
              parameters = list(cutOffs = c(0,1,0), scoring = c(0.25,0.5,2)),
              parallel = NULL, # parallel = N with N number of cores to use or
              ## a list with cores in the first and machines in the second entry
              ## like list(cores=c(2,4,8), machines=c("bionform1", "bioinform2",
              ## "bioinform3"))
              parallel2 = 1,
              selection = c("t"), # can be "t" or "s"
              relFit = FALSE,
              method = "s",
              type = "SOCK",
              exhaustive = FALSE,
              delcyc = TRUE,
              ...
              ) {
        addPriorKnowledge <- get("addPriorKnowledge",
                                 envir = asNamespace("CellNOptR"))
        method <- checkMethod(method)
        if (parameters$cutOffs[1] > parameters$cutOffs[2]) {
            parameters$cutOffs <- sort(parameters$cutOffs)
            print(paste("your're cutoff parameters didn't make any sense.
I can't let you do this, Dave. I changed them to ", parameters$cutOffs, ".",
sep = ""))
        }
        if (is.null(elitism) == TRUE) { elitism <- ceiling(popSize*0.1) }
        if (elitism >= popSize) { elitism <- floor(0.1*popSize) }
        if (is.null(inversion) == TRUE) { inversion <- 0 }
        if (is.null(initBstring) == TRUE) {
            initBstring <- rep(1, length(model$reacID))
        }
        if (length(initBstring) == 1 & length(model$reacID) > 1) {
            initBstring <- max(0, min(initBstring, floor((popSize)*0.5)))
            initBstring <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist,
                                       model=model, approach=approach,
                                       seeds=initBstring, parameters=parameters,
                                       verbose = verbose, parallel=parallel,
                                       parallel2=parallel2,relFit = relFit,
                                       method = method, ...)
            localString <- initBstring
            initBstring <- initBstring[order(localString[, ncol(localString)],
                                             decreasing = TRUE),
                                       -ncol(initBstring)]
        }
        if (!is(CNOlist, "CNOlist")) {
            CNOlist = CNOlist(CNOlist)
        }
        ## create foldchanges if not already included in nemlist
        NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist,
                                parameters = parameters, approach = approach,
                                method)
        spaceExp <- 2^length(model$reacID)
        if ((popSize*stallGenMax) >= spaceExp & exhaustive) {
            print(paste("the genetic algorithm would score at least ",
                        popSize*stallGenMax,
                        " networks, while the size of
the search space is only ", spaceExp,
"; therefore an exhaustive search is initialised", sep = ""))
            result <- exSearch(CNOlist,model,sizeFac,NAFac,NEMlist,parameters,
                               parallel=parallel,relFit = relFit,
                               method = method, ...)
            PopTolScores <-
                result$scores[which(result$scores <
                                    (result$score + abs(result$score)*relTol))]
            PopTol <-
                result$bStrings[which(result$scores <
                                      (result$score +
                                       abs(result$score)*relTol))]
            return(list(bString = result$bString, stringsTol = PopTol,
                        stringsTolScores = PopTolScores,
                        dtmRatio = result$dtmRatio, population = result$scores))
        } else {
            if ("s" %in% selection) {
                if (length(selPress) == 2) {
                    selPressPct <- selPress[2]
                    selPress <- selPress[1]
                } else {
                    selPressPct <- 0
                }
                if (selPress < 1) {
                    print("with selPress less than 1
low ranking networks are favoured")
                    msg1 <- 1
                }
                if (selPress > 2 & fit == "linear") {
                    selPress <- 2
                    print("if selPress is greater than 2,
fit cannot be set to linear; selPress set to 2; or restart with
fit set to nonlinear")
                    msg2 <- 1
                }
                if (selPress >= popSize) {
                    selPress <- popSize - 1
                    print(paste("selPress must be lower than popSize;
selPress set to ", selPress, sep = ""))
                    msg3 <- 1
                }
                if (selPress < 0) {
                    selPress <- 1
                    print("selPress has to be greater than zero and should be
greater than 1; selPress set to 1")
                    msg3 <- 1
                }
            }
            if ("t" %in% selection) {
                selPress <- floor(min(max(2, selPress[1]), popSize/2))
            }
            bLength <- length(model$reacID) # changed from length(initBstring)
            ## simList = prep4sim(model)
            simList <- NULL
            indexList = indexFinder(CNOlist, model)
            ## initialize starting population
            if (is.null(nrow(initBstring))) {
                initBstring <- t(as.matrix(initBstring))
            }
            Pop <- rbind(1 - initBstring, matrix(sample(c(0,1),
            (bLength * (popSize - (nrow(initBstring)*2))), replace = TRUE),
            nrow = (popSize - (nrow(initBstring)*2)), ncol = bLength),
            initBstring) # also the inverted initBstring is added
            Pop <- addPriorKnowledge(Pop, priorBitString)
            bestbit <- Pop[1, ]
            bestobj <- Inf
            stop <- FALSE
            g <- 0
            stallGen <- 0
            res <- rbind(c(g, bestobj, toString(bestbit), stallGen, Inf,
                           Inf, toString(bestbit), 0),
                         c(g, bestobj, toString(bestbit),
                           stallGen, Inf, Inf, toString(bestbit), 0))
            colnames(res) <- c("Generation", "Best_score", "Best_bitString",
                               "Stall_Generation", "Avg_Score_Gen",
                               "Best_score_Gen",
                               "Best_bit_Gen", "Iter_time")
            ## introduce the outputvector based on elitism on or off
            if (elitism >= 1) {
                res <- rbind(c(paste(res[4], " / ", res[1], sep = ""),
                               paste(res[6], " (", res[5], ")", sep = ""),
                               res[7], res[8], "."))
                colnames(res) <- c("Stall Generations / Total Generations",
                                   "Best Score (Average Score)",
                                   "Best_bitString", "Iter_time",
                                   "----------------------------------")
            } else {
                res <- rbind(c(paste(res[4], " / ", res[1], sep = ""),
                               paste(res[6], " (", res[5], ")", sep = ""),
                               res[7], res[2], res[3],res[8], "."))
                colnames(res) <- c("Stall Generations / Total Generations",
                                   "Best Score (Average Score)_Gen",
                                   "Best_bitString_Gen", "Best Score",
                                   "Best_bitString","Iter_time",
                                   "----------------------------------")
            }
            PopTol <- rep(NA, bLength)
            PopTolScores <- NA
            getObj <- function(x) {
                Score <- computeScoreNemT1(CNOlist=CNOlist, model=model,
                                           bString=x, sizeFac=sizeFac,
                                           NAFac=NAFac, approach = approach,
                                           NEMlist = NEMlist,
                                           parameters=parameters, tellme = 0,
                                           relFit = relFit, method = method,
                                           ...)
                return(Score)
            }
            t0 <- Sys.time()
            t <- t0
            selPressOrg <- selPress
            ## add a graph that shows the improvement:
            if (graph) {
                distances <- numeric()
                for (i in seq_len(100)) {
                    distances <- c(distances, dist(rbind(rep(1, i), rep(0, i))))
                }
                graphVal <- numeric()
                graphValAvg <- numeric()
                graphValSp <- numeric()
                gUp <- 0
                graphAxis <- numeric()
                ## graphics.off()
                par(bg = "white")
                graphValDist <- numeric()
                bestMemTurn <- 0
                bestMem <- "off"
            }
            if (!is.null(parallel)) {
                if (is.list(parallel)) {
                    if (length(parallel[[1]]) != length(parallel[[2]])) {
                        stop("The nodes (second list object in parallel) and
the number of cores used on every node (first list object in parallel) must be
the same.")
                    }
                    hosts <- character()
                    for (i in seq_len(length(parallel[[1]]))) {
                        hosts <- c(hosts, rep(parallel[[2]][i],
                                              parallel[[1]][i]))
                    }
                    hosts <- as.list(hosts)
                    sfInit(parallel=TRUE, socketHosts=hosts)
                } else {
                    sfInit(parallel=TRUE, cpus=parallel, type = type)
                }
                sfLibrary(bnem)
                ## sfExport(list = exportVars("ga"), local = TRUE)
            }
            scores <- numeric(nrow(Pop))
            if ("t" %in% selection) {
                tsReduce <- function(x, scores) {
                    y <- x[order(scores[x])[1]]
                    return(y)
                }
            }
            if ("s" %in% selection) {
                susSel <- function(x, wheel1, breaks) {
                    y <- which(wheel1 > breaks[x])[1]
                    return(y)
                }
            }
            popMem <- numeric() # see high diversity later
            scoresMem <- numeric() # see high diversity later
            likelihoods <- numeric()
            lh.samples <- character()
            while (!stop) {
                if (delcyc) {
                    if (!is.null(parallel)) {
                        Pop <- t(sfApply(Pop, 1, removeCycles, model))
                    } else {
                        Pop <- t(apply(Pop, 1, removeCycles, model))
                    }
                }
                bestReduced <- reduceGraph(Pop[popSize, ], model, CNOlist)
                if (sum((Pop[popSize, ] - bestReduced) != 0) > 0) {
                    Pop[1, ] <- bestReduced
                }
                scores <- numeric(nrow(Pop))
                if (!is.null(parallel)) {
                    scores <- sfApply(Pop, 1, getObj)
                } else {
                    scores <- apply(Pop, 1, getObj)
                }
                scores[which(is.nan(scores) == TRUE)] <-
                    max(scores[which(is.nan(scores) == FALSE)]) + 1
                rankP <- order(scores, decreasing = TRUE)
                Pop <- Pop[rankP, ]
                ## something to remember the samples:
                ## for (sample in 1:nrow(Pop)) {
                ##   if (!(toString(Pop[sample, ]) %in% lh.samples)) {
                ##     likelihoods <- c(likelihoods, scores[sample])
                ##     lh.samples <- c(lh.samples, toString(Pop[sample, ]))
                ##   }
                ## }
                ## try to alternatively keep in mind the best networks:
                if (graph) {
                    if (bestMem == "off") {
                        bestMem <- list()
                        bestMem[[1]] <- Pop[popSize, ]
                        bestMem[[2]] <- Pop[popSize, ]
                        bestMemTurn <- abs(bestMemTurn - 1)
                    } else {
                        bestMem[[(bestMemTurn+1)]] <- Pop[popSize, ]
                        bestMemTurn <- abs(bestMemTurn - 1)
                    }
                }
                ## replace worst ranked networks with random networks for
                ## variation, ranking again would slow up the process and
                ## the idea is not to get random good networks but just good
                ## chunks in overall worse networks
                if (inversion >= 1) {
                    Pop[seq_len(inversion), ] <-
                        abs(Pop[(popSize - inversion + 1):popSize, ] - 1)
                }
                scores <- scores[rankP]
                ## for high diversity, similar to sokolov & whitley (2005)
                ## popFull <- rbind(popMem, Pop)
                ## scoresFull <- c(scoresMem, scores)
                ## rankPFull <- order(scoresFull, decreasing = TRUE)
                ## popFull <- popFull[rankPFull, ]
                ## scoresFull <- scoresFull[rankPFull]
                ## bitsInDec <- apply(popFull, 1, bitToDec)
                ## uniquePop <- unique(bitsInDec)
                ## uniquePopPos <- numeric()
                ## for (i in length(scoresFull):(1 +
                ## max(popSize,(length(scoresFull) - length(uniquePop))))) {
                ##   uniquePopPos <- c(uniquePopPos,
                ## which(bitsInDec %in% uniquePop[i])[1])
                ## }
                ## uniquePopPos <- c(uniquePopPos, sample(1:nrow(popFull),
                ## (popSize - length(uniquePopPos))))
                ## popMem <- popFull[uniquePopPos, ]
                ## scoresMem <- scores[uniquePopPos]
                ## Pop <- popMem
                ## selection (s or t)
                PSize3 <- popSize - elitism
                if ("s" %in% selection) {
                    ## nonlinear or linear fitness
                    if (fit == "nonlinear") {
                        x <- as.double(polyroot(c(rep(selPress, popSize-1),
                        (selPress - popSize))))
                        y <- polyroot(c(rep(selPress, popSize-1),
                        (selPress - popSize)))
                        x <- x[order(abs(y - as.double(y)))[1]]
                        fitness <-
                            (popSize*x^(seq_len(popSize-1)))/(
                                sum(x^(seq_len(popSize-1))))
                    }
                    if (fit == "linear") {
                        fitness <- 2 - selPress + (2 * (selPress - 1) *
                                                   (c(seq_len(popSize))
                                                       - 1)/(popSize - 1))
                    }
                    wheel1 <- cumsum(fitness/sum(fitness))
                    breaks <- runif(1) * 1/popSize
                    breaks <- c(breaks, breaks +
                                        ((seq_len((popSize - 1))))/popSize)
                    sel <- rep(1, popSize)

                    if (!is.null(parallel) & popSize > 10000) {
                        sel <- sfApply(as.matrix(seq_len(popSize)), 1, susSel,
                                       wheel1, breaks)
                    } else {
                        sel <- apply(as.matrix(seq_len(popSize)), 1, susSel,
                                     wheel1, breaks)
                    }
                }
                if ("t" %in% selection) {

                    pRanks <- sample(seq_len(popSize), popSize)
                    t.size <- min(popSize/2, selPress)
                    ppRanks <- matrix(0, popSize, t.size)
                    for (i in seq_len(t.size)) {
                        if (i == 1) {
                            ppRanks[, i] <- pRanks[c(popSize,
                                                     seq_len((popSize-1)))]
                        } else {
                            ppRanks[, i] <- ppRanks[c(popSize,
                                                      seq_len((popSize-1))),
                            (i-1)]
                        }
                    }

                    if (!is.null(parallel) & popSize > 10000) {
                        sel <-
                            as.vector(unlist(sfApply(
                                cbind(pRanks, ppRanks), 1, tsReduce, scores)))
                    } else {
                        sel <- as.vector(unlist(apply(
                            cbind(pRanks, ppRanks), 1, tsReduce, scores)))
                    }

                }
                if ("r" %in% selection) {

                    sel <- rep(seq_len(popSize),
                               round((seq_len(popSize))/sum(seq_len(popSize))*
                                     popSize))

                    if (length(sel) < popSize) {
                        sel <- c(sel, sel[length(sel)])
                    }

                }
                if ("f" %in% selection) {

                    scoresF <- scores*(-1)
                    scoresF <- scoresF - min(scoresF)

                    sel <- rep(seq_len(popSize),
                               round((scoresF)/sum(seq_len(scoresF))*popSize))

                    if (length(sel) < popSize) {
                        sel <- c(sel, sel[seq_len((popSize-length(sel)))])
                    }

                }
                ##print(sel)
                ##print(scores)
                Pop2 <- Pop[sel, ]
                PSize2 <- dim(Pop2)[1]
                mates <- cbind(ceiling(runif(PSize3) * PSize2),
                               ceiling(runif(PSize3) * PSize2))
                InhBit <- matrix(runif((PSize3 * bLength)), nrow = PSize3,
                                 ncol = bLength)
                InhBit <- InhBit < 0.5
                Pop3par1 <- Pop2[mates[, 1], ]
                Pop3par2 <- Pop2[mates[, 2], ]
                Pop3 <- Pop3par2
                Pop3[InhBit] <- Pop3par1[InhBit]
                MutProba <- matrix(runif((PSize3 * bLength)), nrow = PSize3,
                                   ncol = bLength)
                MutProba <- (MutProba < (pMutation/bLength))
                Pop3[MutProba] <- 1 - Pop3[MutProba]
                t <- c(t, Sys.time())
                g <- g + 1
                thisGenBest <- scores[popSize]
                thisGenBestBit <- Pop[popSize, ]
                if (is.na(thisGenBest)) {
                    thisGenBest <- min(scores, na.rm = TRUE)
                    thisGenBestBit <- Pop[which(scores == thisGenBest)[1],
                                          ]
                }
                if (thisGenBest < bestobj) {
                    bestobj <- thisGenBest
                    bestbit <- thisGenBestBit
                    stallGen <- 0
                }
                else {
                    stallGen <- stallGen + 1
                }
                if ("s" %in% selection) {
                    if (selPressPct < 0) {
                        selPress <- max(selPress - selPress*selPressPct, 1)
                    } else {
                        if (fit == "linear") {
                            selPress <- min(selPress + selPress*selPressPct, 2)
                        } else {
                            selPress <- min(selPress + selPress*selPressPct,
                                            popSize-1)
                        }
                    }
                }
                ## the following code needs changes if elitism is set to 0:
                resThisGen <- c(g, bestobj, toString(bestbit), stallGen,
                (mean(scores, na.rm = TRUE)), thisGenBest,
                toString(thisGenBestBit),
                as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"),
                "----------------------------------")
                names(resThisGen) <- c("Generation", "Best_score",
                                       "Best_bitString",
                                       "Stall_Generation", "Avg_Score_Gen",
                                       "Best_score_Gen",
                                       "Best_bit_Gen", "Iter_time",
                                       "----------------------------------")
                ## verbose output based on elitism on or off
                if (elitism >= 1) {
                    resThisGen <- c(g, bestobj, toString(bestbit), stallGen,
                    (mean(scores, na.rm = TRUE)), thisGenBest,
                    toString(thisGenBestBit),
                    as.numeric((t[length(t)] - t[length(t) - 1]),
                               units = "secs"),
                    "----------------------------------")
                    names(resThisGen) <- c("Generation", "Best_score",
                                           "Best_bitString",
                                           "Stall_Generation", "Avg_Score_Gen",
                                           "Best_score_Gen",
                                           "Best_bit_Gen", "Iter_time",
                                           "----------------------------------")
                    if (targetBstring[1] == "none") {
                        resThisGen <- c(paste(resThisGen[4], " / ",
                                              resThisGen[1], sep = ""),
                                        paste(resThisGen[6], " (",
                                              resThisGen[5], ")", sep = ""),
                                        resThisGen[7], resThisGen[8],
                                        "----------------------------------")
                        names(resThisGen) <-
                            c("Stall Generations / Total Generations",
                              "Best Score (Average Score)", "Best_bitString",
                              "Iter_time", "----------------------------------")
                    } else {
                        resThisGen <- c(paste(resThisGen[4], " / ",
                                              resThisGen[1], sep = ""),
                                        paste(resThisGen[6], " (",
                                              resThisGen[5], ")", sep = ""),
                                        resThisGen[7], paste(targetBstring,
                                                             collapse = ", "),
                                        resThisGen[8],
                                        "----------------------------------")
                        names(resThisGen) <-
                            c("Stall Generations / Total Generations",
                              "Best Score (Average Score)", "Best_bitString",
                              "Target_bitString", "Iter_time",
                              "----------------------------------")
                    }
                } else {
                    resThisGen <- c(g, bestobj, toString(bestbit), stallGen,
                    (mean(scores, na.rm = TRUE)), thisGenBest,
                    toString(thisGenBestBit),
                    as.numeric((t[length(t)] - t[length(t) - 1]),
                               units = "secs"),
                    "----------------------------------")
                    names(resThisGen) <- c("Generation", "Best_score",
                                           "Best_bitString",
                                           "Stall_Generation",
                                           "Avg_Score_Gen", "Best_score_Gen",
                                           "Best_bit_Gen", "Iter_time",
                                           "----------------------------------")
                    if (targetBstring[1] == "none") {
                        resThisGen <- c(paste(resThisGen[4], " / ",
                                              resThisGen[1], sep = ""),
                                        paste(resThisGen[6], " (",
                                              resThisGen[5], ")", sep = ""),
                                        resThisGen[7], resThisGen[2],
                                        resThisGen[3],resThisGen[8],
                                        "----------------------------------")
                        names(resThisGen) <-
                            c("Stall Generations / Total Generations",
                              "Best Score (Average Score)_Gen",
                              "Best_bitString_Gen", "Best Score",
                              "Best_bitString","Iter_time",
                              "----------------------------------")
                    } else {
                        resThisGen <- c(paste(resThisGen[4], " / ",
                                              resThisGen[1], sep = ""),
                                        paste(resThisGen[6], " (",
                                              resThisGen[5], ")", sep = ""),
                                        resThisGen[7], resThisGen[2],
                                        resThisGen[3], paste(targetBstring,
                                                             collapse = ", "),
                                        resThisGen[8],
                                        "----------------------------------")
                        names(resThisGen) <-
                            c("Stall Generations / Total Generations",
                              "Best Score (Average Score)_Gen",
                              "Best_bitString_Gen", "Best Score",
                              "Best_bitString", "Target_bitString",
                              "Iter_time", "----------------------------------")
                    }
                }
                if (verbose) {
                    if (elitism >= 1) {
                        if (targetBstring[1] == "none") {
                            print(paste(names(resThisGen)[3], ":", sep = ""))
                            print(as.vector(resThisGen[3]))
                            print(paste(names(resThisGen)[1], ": ",
                                        resThisGen[1], sep = ""))
                            print(paste(names(resThisGen)[2], ": ",
                                        resThisGen[2], sep = ""))
                            print(paste(names(resThisGen)[4], ": ",
                                        paste(resThisGen[4], "s @ ",
                                              Sys.time(), sep = ""), sep = ""))
                            print(paste(names(resThisGen)[5], ": ",
                                        resThisGen[5], sep = ""))
                        } else {
                            print(paste(names(resThisGen)[3], ":", sep = ""))
                            print(as.vector(resThisGen[3]))
                            print(paste(names(resThisGen)[4], ":", sep = ""))
                            print(as.vector(resThisGen[4]))
                            print(paste(names(resThisGen)[1], ": ",
                                        resThisGen[1], sep = ""))
                            print(paste(names(resThisGen)[2], ": ",
                                        resThisGen[2], sep = ""))
                            print(paste(names(resThisGen)[5], ": ",
                                        resThisGen[5], sep = ""))
                            print(paste(names(resThisGen)[6], ": ",
                                        paste(resThisGen[6], "s @ ",
                                              Sys.time(), sep = ""), sep = ""))
                        }
                    } else {
                        print(resThisGen)
                    }
                }
                res <- rbind(res, resThisGen)
                Criteria <- c((stallGen > stallGenMax),
                (as.numeric((t[length(t)] -
                             t[1]), units = "secs") > maxTime), (g > maxGens))
                ## introduce a stop criteria for a target network
                if (targetBstring[1] != "none" & g >= 2) {
                    Criteria <- c((stallGen > stallGenMax),
                    (as.numeric((t[length(t)] -
                                 t[1]), units = "secs") > maxTime),
                    (g > maxGens),
                    (computeScoreNemT1(CNOlist=CNOlist,
                                       model=model,
                                       bString=bestbit,
                                       sizeFac=sizeFac,
                                       NAFac=NAFac,
                                       approach = approach,
                                       NEMlist = NEMlist,
                                       parameters=parameters, tellme = 0,
                                       relFit = relFit, ...) <=
                     computeScoreNemT1(CNOlist=CNOlist, model=model,
                                       bString=targetBstring, sizeFac=sizeFac,
                                       NAFac=NAFac, approach = approach,
                                       NEMlist = NEMlist, parameters=parameters,
                                       tellme = 0, relFit = relFit, ...)))
                }
                if (any(Criteria)) {
                    stop <- TRUE
                }
                tolScore <- abs(scores[popSize] * relTol)
                TolBs <- which(scores < scores[popSize] + tolScore)
                if (length(TolBs) > 0) {
                    PopTol <- rbind(PopTol, Pop[TolBs, ])
                    PopTolScores <- c(PopTolScores, scores[TolBs])
                }
                if (elitism > 0) {
                    Pop <- rbind(Pop3, Pop[(popSize - elitism + 1):popSize, ])
                }
                else {
                    Pop <- Pop3
                }
                Pop <- addPriorKnowledge(Pop, priorBitString)
                if (graph) {
                    split.screen(figs = c(1, 2), erase = FALSE)
                    split.screen(figs = c(2, 1), screen = 1, erase = FALSE)
                    graphDraw <- 0
                    if (length(graphVal) > 0) {
                        if ((scores[order(scores, decreasing = TRUE)[popSize]] <
                             graphVal[length(graphVal)]) |
                            (stallGen %in% ceiling(stallGenMax*c(1)))) {
                            ModelCut <- model
                            ModelCut$interMat <-
                                ModelCut$interMat[, as.logical(Pop[popSize, ])]
                            ModelCut$notMat <-
                                ModelCut$notMat[, as.logical(Pop[popSize, ])]
                            ModelCut$reacID <-
                                ModelCut$reacID[as.logical(Pop[popSize, ])]
                            screen(2, new = TRUE)
                            plotDnf(ModelCut$reacID, CNOlist = CNOlist)
                                        #}
                            graphVal <-
                                c(graphVal,
                                  scores[order(scores,
                                               decreasing = TRUE)[popSize]])
                            graphValAvg <- c(graphValAvg,
                                             sum(scores)/length(scores))
                            graphValSp <- c(graphValSp, selPress)
                            gUp <- gUp + 1
                            graphAxis <- c(graphAxis, (g+1))
                            graphDraw <- 1
                            graphValDist <- c(graphValDist,
                                              dist(rbind(bestMem[[1]],
                                                         bestMem[[2]])))
                        }
                    } else {
                        graphVal <-
                            c(graphVal, scores[order(scores,
                                                     decreasing =
                                                         TRUE)[popSize]])
                        graphValAvg <-
                            c(graphValAvg, sum(scores)/length(scores))
                        graphValSp <- c(graphValSp, selPress)
                        graphValDist <-
                            c(graphValDist, dist(rbind(bestMem[[1]],
                                                       bestMem[[2]])))
                        gUp <- gUp + 1
                        graphAxis <- c(graphAxis, (g+1))
                        graphDraw <- 1
                    }
                    if (graphDraw == 1) {
                        screen(3, new = TRUE)
                        plot(seq_len(gUp), graphVal, col = "red", type = "l",
                             main = paste("Score Improvement", sep = ""),
                             xlab = "Generation", ylab = "Score",
                             ylim = c(min(graphVal),
                                      max(c(graphVal, graphValAvg))),
                             xaxt = "n")
                        axis(1, at = seq_len(gUp), labels = graphAxis)
                        lines(seq_len(gUp), graphValAvg, col = "blue",
                              type = "l")
                        if (selection %in% "s") {
                            legend(gUp, max(c(graphVal, graphValAvg)),
                                   legend = c("Best Score", "Average Score",
                                              "Selective Pressure"),
                                   fill = c("red", "blue", "black"), xjust = 1)
                            mtext("Selective Pressure", 4, line = 2)
                            par(new=TRUE)
                            plot(seq_len(gUp), graphValSp, col = "black",
                                 type = "l", axes=FALSE, ylab = "", xlab = "",
                                 yaxt = "n",
                                 ylim = c(1,max(2,ceiling(max(graphValSp)))))
                            axis(4, at = (10:max(20,(ceiling(max(
                                                 graphValSp))*10)))/10,
                                 labels = (10:max(20,(ceiling(max(
                                                  graphValSp))*10)))/10)
                        } else {
                            legend(gUp, max(c(graphVal, graphValAvg)),
                                   legend = c("Best Score", "Average Score"),
                                   fill = c("red", "blue"), xjust = 1)
                        }
                        screen(4, new = TRUE)
                        plot(seq_len(gUp), graphValDist, col = "black",
                             type = "l", main = "Best Strings Distance",
                             xlab = "Generation", ylab = "Distance",
                             ylim = c(min(graphValDist), max(graphValDist)),
                             xaxt='n')
                        mtext("Edge Difference", 4, line = 2)
                        abline(h = distances, lty = 3)
                        axis(4, at = distances,
                             labels = seq_len(length(distances)))
                        axis(1, at = seq_len(gUp), labels = graphAxis)
                    }
                }
            }
            PopTol <- as.matrix(PopTol)[-1, ]
            PopTolScores <- PopTolScores[-1]
            TolBs <- which(PopTolScores < scores[popSize] + tolScore)
            PopTol <- as.matrix(PopTol)[TolBs, ]
            PopTolScores <- PopTolScores[TolBs]
            PopTolT <- cbind(PopTol, PopTolScores)
            PopTolT <- unique(PopTolT, MARGIN = 1)
            if (!is.null(dim(PopTolT))) {
                PopTol <- PopTolT[, seq_len((dim(PopTolT)[2] - 1))]
                PopTolScores <- PopTolT[, dim(PopTolT)[2]]
            }
            else {
                PopTol <- PopTolT[seq_len((length(PopTolT) - 1))]
                PopTolScores <- PopTolT[length(PopTolT)]
            }
            res <- res[2:dim(res)[1], ]
            rownames(res) <- NULL
            if (!is.null(parallel)) {
                sfStop()
            }
            bestbit <- reduceGraph(bestbit, model, CNOlist)
            dtmRatio <- 0
            return(list(bString = bestbit, results = res, stringsTol = PopTol,
                        stringsTolScores = PopTolScores, dtmRatio = dtmRatio,
                        population = likelihoods))
        }
    }
#' @noRd
#' @importFrom mnem plotDnf
getHierarchy <-
    function(graph) {
        adj <- dnf2adj(graph)
        dnf <- adj2dnf(adj)
        g <- plotDnf(dnf, draw = FALSE)
        Ypos <- g@renderInfo@nodes$labelY
        Ynames <- names(g@renderInfo@nodes$labelY)
        ## Ypos <- Ypos[-grep("and", Ynames)]
        ## Ynames <- Ynames[-grep("and", Ynames)]
        hierarchy <- list()
        count <- 0
        for (i in sort(unique(Ypos), decreasing = TRUE)) {
            count <- count + 1
            hierarchy[[count]] <- Ynames[which(Ypos == i)]
        }
        return(hierarchy)
    }
#' @noRd
#' @import
#' matrixStats
#' stats
#' @importFrom flexclust dist2
getNemFit <-
    function (simResults, CNOlist, model, indexList, # VERSION of CNO: 1.4
              timePoint = c("t1", "t2"),
              NAFac = 1, sizePen, simResultsT0 = NA,
              tellme = 0,
              approach = "fc",
              NEMlist,
              parameters,
              sim = 0,
              relFit = FALSE,
              method = "pearson",
              verbose = FALSE,
              opt = "min"
              ) {
        if (!is(CNOlist, "CNOlist")) {
            CNOlist = CNOlist(CNOlist)
        }
        simResults <- simResults[, colnames(CNOlist@signals[[1]])]
        if (timePoint == "t1") {
            tPt <- 2
        }
        else {
            if (timePoint == "t2") {
                tPt <- 3
            }
            else {
                tPt <- timePoint
            }
        }
        ## MY CODE START start with the nem approaches for the error
        MSEAabs <- NULL
        MSEIabs <- NULL
        MSEAfc <- NULL
        MSEIfc <- NULL
        NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist,
                                parameters = parameters, approach = approach,
                                method = method)
        if ("abs" %in% approach) {
            MSEE <- rep(Inf, nrow(NEMlist$norm))
            S.mat <- simResults
            colnames(NEMlist$norm)[which(colnames(NEMlist$norm) %in% "")] <-
                "Base"
            rownames(S.mat)[which(rownames(S.mat) %in% "")] <- "Base"
            S.mat <- S.mat[colnames(NEMlist$norm), ]
            if (any(c("spearman", "pearson", "kendall") %in% method)) {
                if ("spearman" %in% method) {
                    E.mat <- NEMlist$exprs
                    S.mat.ranks <- t(S.mat)
                    cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p"))
                } else {
                    E.mat <- NEMlist$exprs
                    cosine.sim <- -t(cor(S.mat, t(E.mat), method = method))
                }
                R <- cbind(cosine.sim, -cosine.sim)
            } else {
                ESpos <- NEMlist$norm%*%abs(1 - S.mat)*parameters$scoring[1]
                ES0 <- abs(1 - NEMlist$norm)%*%S.mat*parameters$scoring[2]
                ESposI <- NEMlist$norm%*%S.mat*parameters$scoring[1]
                ES0I <- abs(1 -
                            NEMlist$norm)%*%abs(1 - S.mat)*parameters$scoring[2]
                MSEAabs <- (ESpos + ES0)
                MSEIabs <- (ESposI + ES0I)
                R <- cbind(MSEAabs, MSEIabs)
                R <- R/ncol(NEMlist$exprs)
            }
            R[is.na(R)] <- max(R[!is.na(R)])
            MSEE <- rowMins(R)
        }
        if ("fc" %in% approach) {
            MSEE <- rep(Inf, nrow(NEMlist$fc))
            SCompMat <- computeFc(CNOlist, t(simResults))
            signtmp <- sign(SCompMat)
            SCompMat <- SCompMat/SCompMat
            SCompMat[is.na(SCompMat)] <- 0
            SCompMat <- SCompMat*signtmp
            SCompMat <- t(SCompMat)
            ## for debugging:
            ## if (any(!(colnames(NEMlist$fc) %in% rownames(SCompMat)))) {
            ##     print(colnames(NEMlist$fc)[which(!(colnames(NEMlist$fc) %in%
            ##                                        rownames(SCompMat)))]);
            ## }
            if (is.null(rownames(SCompMat))) {
                SCompMat <- SCompMat[, colnames(NEMlist$fc)]
            } else {
                SCompMat <- SCompMat[colnames(NEMlist$fc), ]
            }
            if (ncol(SCompMat) != ncol(NEMlist$fc)) {
                SCompMat <- t(SCompMat)
            }
            if (length(grep("_vs_", rownames(SCompMat))) > 0) {
                SCompMat <- t(SCompMat)
            }
            if (any(c("cosine", "euclidean", "maximum", "manhattan", "canberra",
                      "binary", "minkowski", "spearman", "pearson",
                      "kendall") %in% method)) {
                S.mat <- SCompMat
                E.mat <- NEMlist$fc
                if (any(c("euclidean", "maximum", "manhattan", "canberra",
                          "binary", "minkowski") %in% method)) {
                    power <- as.numeric(method)
                    power <- power[-which(is.na(power)==TRUE)]
                    if (length(power) == 0) {
                        power <- 2
                    }
                    method <- method[which(method %in% c("euclidean", "maximum",
                                                         "manhattan",
                                                         "canberra", "binary",
                                                         "minkowski"))[1]]
                    cosine.sim <-  dist2(S.mat, E.mat, method = method,
                                         p = power)
                    cosine.simI <- dist2(S.mat, -E.mat, method = method,
                                         p = power)
                    MSEAfc <- t(cosine.sim)
                    MSEIfc <- t(cosine.simI)
                }
                if (any(c("spearman", "pearson", "kendall",
                          "cosine") %in% method)) {
                    weighted <- parameters$scoring[2]
                    if (weighted > 1) {
                        for (i in seq_len(nrow(S.mat))) {
                            S.mat[i, grep(rownames(S.mat)[i],
                                          colnames(S.mat))] <-
                                S.mat[i, grep(rownames(S.mat)[i],
                                              colnames(S.mat))]*weighted
                        }
                    }
                    if (parameters$scoring[1] == 0) {
                        S.mat[which(S.mat == 0)] <- NA
                    }
                    if ("pearson" %in% method) {
                        cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p",
                                             use = "pairwise.complete.obs"))
                    }
                    if ("spearman" %in% method) {
                        S.mat.ranks <- t(S.mat)
                        cosine.sim <- -t(cor(S.mat.ranks, t(E.mat),
                                             method = "p",
                                             use = "pairwise.complete.obs"))
                    }
                    if ("kendall" %in% method) {
                        cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "k",
                                             use = "pairwise.complete.obs"))
                    }
                    if ("cosine" %in% method) {
                        if (any(is.na(S.mat))) {
                            S.mat[is.na(S.mat)] <- 0
                        }
                        vprod <- function(x) { return(x%*%x) }
                        S.mat <- S.mat/(apply(S.mat, 1, vprod)^0.5)
                        E.mat <- E.mat/(apply(E.mat, 1, vprod)^0.5)
                        cosine.sim <- -E.mat%*%t(S.mat)
                    }
                    if ("test" %in% method) {
                        cs.fisher <- atanh(cosine.sim)
                        ## 0.5*log((1+cosine.sim)/(1-cosine.sim))
                        cs.sd <- 1/sqrt(ncol(NEMlist$fc) - 3)
                        cs.z <- (cs.fisher - 0)/cs.sd
                        cs.sign <- sign(cosine.sim)
                        cosine.sim <- 2*pnorm(-abs(cs.z))
                        cosine.sim[seq_len(length(cosine.sim))] <-
                            abs(1 - p.adjust(cosine.sim, method = "bonferroni"))
                        cosine.sim <- cosine.sim*cs.sign
                    }
                    cosine.sim[is.na(cosine.sim)] <- 0
                    MSEAfc <- cosine.sim
                    MSEIfc <- -cosine.sim
                }
            } else {
                if ("llr" %in% method) {
                    S.mat <- SCompMat
                    E.mat <- NEMlist$fc
                    cosine.sim <- -E.mat%*%t(S.mat)
                    cosine.sim[is.na(cosine.sim)] <- 0
                    MSEAfc <- cosine.sim
                    MSEIfc <- -cosine.sim
                } else {
                    S0 <- as.matrix(1 - abs(SCompMat))
                    Spos <- SCompMat
                    Spos[which(Spos == -1)] <- 0
                    Sneg <- SCompMat
                    Sneg[which(Sneg == 1)] <- 0
                    if ("mLL" %in% method | "cp" %in% method) {
                        Sneg <- abs(Sneg)
                        E0 <- NEMlist$E0
                        Epos <- NEMlist$Epos
                        Eneg <- NEMlist$Eneg
                        EposI <- NEMlist$EposI
                        EnegI <- NEMlist$EnegI
                        ES0 <- E0%*%S0
                        ES0pos <- E0%*%Spos
                        ES0neg <- E0%*%Sneg
                        ESpos <- Epos%*%Spos
                        ESpos0 <- Epos%*%S0
                        ESposneg <- Epos%*%Sneg
                        ESneg <- Eneg%*%Sneg
                        ESneg0 <- Eneg%*%S0
                        ESnegpos <- Eneg%*%Spos
                        ESposI <- EposI%*%Spos
                        ESposI0 <- EposI%*%S0
                        ESposIneg <- EposI%*%Sneg
                        ESnegI <- EnegI%*%Sneg
                        ESnegI0 <- EnegI%*%S0
                        ESnegIpos <- EnegI%*%Spos
                        alpha1 <- parameters$scoring[1]
                        beta2 <- parameters$scoring[2]
                        beta4 <- beta2/3
                        if ("cont" %in% method) {
                            MSEAfc <- exp(ESpos + ESneg)
                            MSEIfc <- exp(ESposneg + ESnegpos)
                        } else {
                            alpha2 <- alpha1/2
                            beta1 <- parameters$scoring[3]
                            beta3 <- beta2/2
                            Alpha <- 1 - alpha1 - alpha2
                            Beta <- 1 - beta1 - beta2 - beta3 - beta4
                            MSEAfc <- Alpha^ES0 * Beta^(ESpos + ESneg) *
                                beta1^(ESposI + ESnegI) *
                                beta2^(ES0pos + ES0neg) *
                                beta3^(ESposIneg + ESnegIpos) *
                                beta4^(ESposneg + ESnegpos) *
                                alpha1^(ESposI0 + ESnegI0) *
                                alpha2^(ESpos0 + ESneg0)
                            MSEIfc <- Alpha^ES0 * beta4^(ESpos + ESneg) *
                                beta3^(ESposI + ESnegI) *
                                beta2^(ES0pos + ES0neg) *
                                beta1^(ESposIneg + ESnegIpos) *
                                Beta^(ESposneg + ESnegpos) *
                                alpha1^(ESposI0 + ESnegI0) *
                                alpha2^(ESpos0 + ESneg0)
                        }
                        if ("cp" %in% method) {
                            MSEAfc <- -log(MSEAfc)
                            MSEIfc <- -log(MSEIfc)
                        }
                    } else {
                        E0 <- NEMlist$E0
                        Epos <- NEMlist$Epos
                        Eneg <- NEMlist$Eneg
                        EposI <- NEMlist$EposI
                        EnegI <- NEMlist$EnegI
                        ## differ between count all fc and alignment
                        ESpos <- Epos%*%Spos
                        ESneg <- Eneg%*%Sneg
                        ESposI <- EposI%*%Spos
                        ESnegI <- EnegI%*%Sneg
                        ES0 <- E0%*%S0
                        MSEAfc <- (ESpos + ESneg + ES0)
                        MSEIfc <- (ESposI + ESnegI + ES0)
                    }
                }
            }
            if ("mLL" %in% method) {
                R <- rbind(MSEAfc, MSEIfc)
                MSEE <- log(rowSums(R))
            } else {
                if (any(c("cosine", "euclidean", "maximum", "manhattan",
                          "canberra", "binary", "minkowski", "spearman",
                          "pearson", "kendall") %in% method)) {
                    R <- cbind(MSEAfc, MSEIfc)
                    R[is.na(R)] <- max(R[!is.na(R)])
                    if (relFit) {
                        S0 <- as.matrix(1 - abs(SCompMat))
                        Spos <- SCompMat
                        Spos[which(Spos == -1)] <- 0
                        Sneg <- SCompMat
                        Sneg[which(Sneg == 1)] <- 0
                        SS0 <- t(S0*parameters$scoring[1])%*%S0
                        SSpos <- t(Spos)%*%Spos
                        SSneg <- t(Sneg)%*%Sneg
                        MSES <- (SSpos + SSneg + SS0)/ncol(NEMlist$fc)
                        R <- t(t(R)*diag(MSES))
                    }
                } else {
                    R <- cbind(MSEAfc, MSEIfc)
                    R[is.na(R)] <- max(R[!is.na(R)])
                    if (relFit) {
                        SS0 <- t(S0*parameters$scoring[1])%*%S0
                        SSpos <- t(Spos)%*%Spos
                        SSneg <- t(Sneg)%*%Sneg
                        MSES <- (SSpos + SSneg + SS0)
                        R <- t(t(R)/diag(MSES))
                        R[is.na(R)] <- 0
                    } else {
                        R <- R#/ncol(NEMlist$fc)
                    }
                }
                if (is.null(NEMlist$egenes)) {
                    MSEE <- rowMins(R)
                } else {
                    MSEE <- numeric()
                    R <- R+NEMlist$geneGrid
                    R[is.na(R)] <- Inf
                    if (is.null(NEMlist$weights)) {
                        MSEE <- rowMins(R)
                    } else {
                        R[is.infinite(R)] <- 0
                        topNsum <- function(x, N) {
                            y <- sum(x[order(x)[seq_len(N)]])
                            return(y)
                        }
                        MSEE <- apply(R, 1, topNsum, ncol(CNOlist@signals[[1]]))
                        MSEE <- MSEE*NEMlist$weights
                    }
                }
            }
        }
        if (parameters$cutOffs[3] == -1) {
            ## do median polish over gene clusters
            data.med <- NEMlist$fc[seq_len(ncol(CNOlist@signals[[2]])), ]*0
            Epos <- which(R == MSEE, arr.ind = TRUE)
            for (i in seq_len(ncol(CNOlist@signals[[1]]))) {
                tmp <-
                    medpolish(rbind(
                        NEMlist$fc[Epos[
                                    which(Epos[, 2] == i), 1], ],
                        -NEMlist$fc[Epos[which(Epos[, 2] ==
                                               (i+ncol(CNOlist@signals[[2]]))),
                                         1], ]), trace.iter=FALSE)
                data.med[i, ] <- tmp$col
            }
            E.mat <- data.med
            E.mat[is.na(E.mat)] <- 0
            tmp <- which(apply(E.mat, 1, sum) != 0)
            E.mat <- as.matrix(E.mat[which(apply(E.mat, 1, sum) != 0), ])
            rownames(E.mat) <- rownames(S.mat)[tmp]
            NEMlist$fc <- E.mat
            S.mat <- SCompMat
            if (parameters$scoring[1] == 0) {
                S.mat[which(S.mat == 0)] <- NA
            }
            if ("pearson" %in% method) {
                cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p",
                                     use = "pairwise.complete.obs"))
            }
            if ("spearman" %in% method) {
                S.mat.ranks <- t(S.mat)
                cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p",
                                     use = "pairwise.complete.obs"))
            }
            cosine.sim[is.na(cosine.sim)] <- 0
            MSEAfc <- cosine.sim
            MSEIfc <- -cosine.sim
            R <- cbind(MSEAfc, MSEIfc)
            R[is.na(R)] <- max(R[!is.na(R)])
            MSEE <- rowMins(R)
        }
        if (tellme == 1) {
            names(MSEE) <- rownames(NEMlist$fc)
            if (parameters$cutOffs[3] > 0 & parameters$cutOffs[3] <= 1) {
                R <- R[which(MSEE < -parameters$cutOffs[3]), ]
                MSEE <- MSEE[which(MSEE < -parameters$cutOffs[3])]
            }
            if (parameters$cutOffs[3] < 0 & parameters$cutOffs[3] > -1) {
                score.quantile <- quantile(MSEE, -parameters$cutOffs[3])
                R <- R[which(MSEE < score.quantile), ]
                MSEE <- MSEE[which(MSEE < score.quantile)]
            }
            if (parameters$cutOffs[3] > 1 & parameters$cutOffs[3] <=
                length(MSEE)) {
                R <- R[order(MSEE)[seq_len(parameters$cutOffs[3])], ]
                MSEE <- MSEE[order(MSEE)[seq_len(parameters$cutOffs[3])]]
            }
            topEgenes <- 0
            subtopo <- matrix(0, nrow = length(MSEE),
                              ncol = ncol(CNOlist@signals[[1]]))
            colnames(subtopo) <- colnames(CNOlist@signals[[1]])
            if (is.null(NEMlist$weights)) {
                Epos <- which(R == MSEE, arr.ind = TRUE)
            } else {
                Epos <-
                    which(NEMlist$geneGrid[
                                    , seq_len(ncol(CNOlist@signals[[1]]))] == 0,
                          arr.ind = TRUE)
            }
            if (length(Epos) != 0) {
                if (is.null(dim(Epos))) {
                    Epos <- t(as.matrix(Epos))
                }
                posReg <- Epos[which(Epos[, 2] <= ncol(CNOlist@signals[[1]])), ]
                if (length(posReg) > 0) {
                    if (length(posReg) > 3) {
                        subtopo[posReg] <- 1
                    } else {
                        subtopo[posReg[1], posReg[2]] <- 1
                    }
                }
                negReg <- Epos[which(Epos[, 2] > ncol(CNOlist@signals[[1]])), ]
                if (length(negReg) > 0) {
                    if (length(negReg) > 3) {
                        negReg[, 2] <- negReg[, 2] - ncol(CNOlist@signals[[1]])
                        subtopo[negReg] <- -1
                    } else {
                        negReg[2] <- negReg[2] - ncol(CNOlist@signals[[1]])
                        subtopo[negReg[1], negReg[2]] <- -1
                    }
                }

                EtoS <- matrix(0, nrow = nrow(Epos), ncol = 4)
                colnames(EtoS) <- c("Egene", "Sgene", "Type", "MSE")
                rownames(EtoS) <- names(MSEE)[Epos[, 1]]

                Epos[which(Epos[, 2] > ncol(subtopo)), 2] <-
                    Epos[which(Epos[, 2] > ncol(subtopo)), 2] - ncol(subtopo)

                EtoS[, 1] <- Epos[, 1]
                EtoS[, 2] <- Epos[, 2]
                EtoS[, 3] <- subtopo[cbind(Epos[, 1], Epos[, 2])]
                EtoS[, 4] <- MSEE[Epos[, 1]]

                EtoS <- EtoS[order(EtoS[, 4], decreasing = FALSE), ]
                sgeneScore <- numeric(ncol(subtopo))
            } else {
                EtoS <- matrix(0, nrow = 1, ncol = 4)
                colnames(EtoS) <- c("Egene", "Sgene", "Type", "MSE")
            }
            if (verbose) {
                for (j in seq_len(ncol(subtopo))) {
                    print(paste(j, ".", colnames(simResults)[j], ": ",
                                sum(EtoS[, 2] == j), sep = ""))
                    print(paste("Activated: ", sum(EtoS[, 2] == j &
                                                   EtoS[, 3] == 1), sep = ""))
                    print(paste("Inhibited: ", sum(EtoS[, 2] == j &
                                                   EtoS[, 3] == -1), sep = ""))
                    print("Summary Score:")
                    print(summary(EtoS[which(EtoS[, 2] == j), 4]))
                }
                dups <- sum(duplicated(rownames(Epos)) == TRUE)
                if (dups > 0) {
                    used <-
                        sum(EtoS[-which(duplicated(rownames(Epos)) == TRUE), 2]
                            %in% seq_len(ncol(subtopo)))
                } else {
                    used <- nrow(EtoS)
                }
                print(paste("Unique genes used: ",
                (used), " (", round((used/nrow(NEMlist$fc))*100, 2), " %)",
                sep = ""))
                print(paste("Duplicated genes: ", dups, sep = ""))
                print("Overall fit:")
                print(summary(EtoS[, 4]))
            }
        }
        if ("mLL" %in% method) {
            deviationPen <- -sum(MSEE)
        } else {
            if (parameters$cutOffs[3] > 1 &
                parameters$cutOffs[3] <= length(MSEE) & tellme == 0) {
                MSEE <- MSEE[order(MSEE, decreasing = FALSE)[
                    seq_len(parameters$cutOffs[3])]]
            }
            if (parameters$cutOffs[3] > 0 & parameters$cutOffs[3] <= 1 &
                tellme == 0) {
                MSEE <- MSEE[which(MSEE < -parameters$cutOffs[3])]
            }
            if (parameters$cutOffs[3] < 0 & parameters$cutOffs[3] > -1 &
                tellme == 0) {
                score.quantile <- quantile(MSEE, -parameters$cutOffs[3])
                MSEE <- MSEE[which(MSEE < score.quantile)]
            }
            if ("cp" %in% method) {
                deviationPen <- sum(MSEE[!is.na(MSEE)])/nrow(NEMlist$fc)
            } else {
                if (parameters$cutOffs[3] > 1 &
                    parameters$cutOffs[3] <= length(MSEE) & tellme == 0) {
                    deviationPen <-
                        1 + sum(MSEE[!is.na(MSEE)])/parameters$cutOffs[3]
                } else {
                    deviationPen <-
                        1 + sum(MSEE[!is.na(MSEE)])#/nrow(NEMlist$fc)
                }
            }
        }
        ## END of my code
        NAPen <- NAFac * length(which(is.na(simResults)))
        if (opt %in% "min") {
            score <- deviationPen + NAPen + sizePen
        } else {
            score <- 1 - deviationPen - NAPen - sizePen
        }
        if (tellme == 1) {
            return(list(EtoS = EtoS, score = score))
        } else {
            return(score)
        }
    }
#' @noRd
#' @import graph
graph2adj <-
    function(gR) {
        adj.matrix <- matrix(0,
                             length(nodes(gR)),
                             length(nodes(gR))
                             )
        rownames(adj.matrix) <- nodes(gR)
        colnames(adj.matrix) <- nodes(gR)
        for (i in seq_len(length(nodes(gR)))) {
            adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
        }

        return(adj.matrix)
    }
#' @noRd
isDag <-
    function(graph = NULL, bString = 0, model = NULL) {
        if (any(bString != 0)) {
            graph <- model$reacID[which(bString == 1)]
        }
        if (!is.null(graph)) {
            adjmat <- dnf2adj(graph)
            ##.order <- apply(adjmat, 1, sum)
            get.order2 <- apply(adjmat, 2, sum)
            adjmat <- adjmat[order(get.order2, decreasing = FALSE),
                             order(get.order2, decreasing = FALSE)]
            if (all(adjmat[lower.tri(adjmat)] == 0)) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            return(TRUE)
        }
    }
#' @noRd
#' @import cluster
kmeansNorm <-
    function(x, k = 2) {
        if (!is.matrix(x)) {
            x <- t(as.matrix(x))
            if (nrow(x) != 1) {
                x <- t(x)
            }
        }
        y <- x
        for (i in seq_len(nrow(x))) {

            if (sd(x[i, ]) == 0) { next() }

            cat('\r', paste(round(i/nrow(x)*100), "%", sep = ""))
            flush.console()

            x.clust <- kmeans(x[i, ], k)
            x.dist <- dist(x[i, ])
            x.sil <- silhouette(x.clust$cluster, x.dist)
            x.norm <- x.sil[, 3]
            if (sum(x[i, which(x.clust$cluster == 1)]) <
                sum(x[i, which(x.clust$cluster == 2)])) {
                x.norm[which(x.clust$cluster == 1)] <-
                    x.norm[which(x.clust$cluster == 1)]*(-1)
            } else {
                x.norm[which(x.clust$cluster == 2)] <-
                    x.norm[which(x.clust$cluster == 2)]*(-1)
            }
            x.norm <- x.norm - min(x.norm)
            x.norm <- x.norm/max(x.norm)
            y[i, ] <- x.norm
        }
        return(y)
    }
#' @noRd
#' @import
#' matrixStats
#' snowfall
#' CellNOptR
#' @importFrom mnem plotDnf
#' @importFrom Biobase rowMin rowMax
localSearch <-
    function(CNOlist, NEMlist, model, approach = "fc", initSeed = NULL,
             seeds = 1,
             parameters = list(cutOffs = c(0,1,0), scoring = c(0.25,0.5,2)),
             sizeFac = 10^-10, NAFac = 1, relTol = 0.01, verbose = TRUE,
             parallel=NULL, parallel2 = 1, relFit = FALSE, method = "s",
             max.steps = Inf, max.time = Inf, node = NULL, absorpII = TRUE,
             draw = TRUE, prior = NULL) {
        if (is.null(prior)) {
            prior <- rep(0, length(model$reacID))
        }
        method <- checkMethod(method)
        debug <- FALSE
        if (verbose %in% "part") {
            verbose2 <- "part"
            verbose <- TRUE
        } else {
            verbose2 <- ""
        }
        if (seeds == "max") {
            seeds <- length(model$reacID)*2
            seeds2 <- "max"
        } else {
            seeds2 <- "none"
        }
        if (!is(CNOlist, "CNOlist")) {
            CNOlist = CNOlist(CNOlist)
        }
        NEMlist <- checkNEMlist(NEMlist, CNOlist, parameters, approach, method)
        bLength <- length(model$reacID)
        ##simList = prep4sim(model)
        indexList = indexFinder(CNOlist, model)
        n = seeds # number of strings to analyse
        bitStrings <- matrix(0, nrow = n, ncol = bLength)
        if (n >= 1) {
            bitStrings[1, ] <- c(0, rep(0, bLength-1))
        }
        if (!is.null(initSeed)) {
            bitStrings[1, ] <- reduceGraph(initSeed, model, CNOlist)
        }
        ## add a few random (but good) strings:
        if (seeds >= 2 & seeds2 != "max") {
            bitStrings[2:seeds, ] <- createCube(ncol(bitStrings), seeds-1)
        }
        bitStringsMem <- numeric()
        bitStringsScores <- numeric()
        if (!is.null(parallel)) {
            if (is.list(parallel)) {
                if (length(parallel[[1]]) != length(parallel[[2]])) {
                    stop("The nodes (second list object in parallel) and the
 number of cores used on every node (first list object in parallel) must be
 the same.") }
                hosts <- character()
                for (i in seq_len(length(parallel[[1]]))) {
                    hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
                }
                hosts <- as.list(hosts)
                sfInit(parallel=TRUE, cpus=sum(parallel[[1]]), type="SOCK",
                       socketHosts=hosts)
            } else {
                sfInit(parallel=TRUE, cpus=parallel)
            }
            sfLibrary(bnem)
            ## sfExport(list = exportVars("loc"), local = TRUE)
        }
        start.time <- Sys.time()
        processSeed <- function(i, bitStrings, CNOlist, model, simList = NULL,
                                indexList, sizeFac, NAFac, approach = approach,
                                NEMlist = NEMlist, parameters, tellme = 0,
                                relFit = relFit, method = method,
                                max.steps = max.steps, node = node, ...) {
            stimuli <- colnames(CNOlist@stimuli)
            inhibitors <- colnames(CNOlist@inhibitors)
            signals <- colnames(CNOlist@signals[[1]])
            step.count <- 0
            row <- i
            new <- FALSE
            bitString <- bitStrings[row, ]
            compMatTmp <- t(t(bitStringsMem) - bitString)
            compScore <- rowMaxs(abs(compMatTmp))
            if (length(compScore) == 0) { compScore <- 1 }
            if (min(compScore) == 0) {
                new <- FALSE
                safeNumber <- 0
                while(!new & safeNumber < 1000) {
                    safeNumber <- safeNumber + 1
                    bitString <- sample(c(0,1), bLength, replace = TRUE)
                    compMatTmp <-
                        t(t(rbind(bitStrings,bitStringsMem)) - bitStringTmp)
                    compScore <- rowMaxs(abs(compMatTmp))
                    if (min(compScore) > 0) {
                        new <- TRUE
                    }
                }
                print(safeNumber)
            }
            fullScore <- computeScoreNemT1(CNOlist=CNOlist, model=model,
                                           bString=bitString, sizeFac=sizeFac,
                                           NAFac=NAFac, approach = approach,
                                           NEMlist = NEMlist,
                                           parameters=parameters, tellme = 0,
                                           relFit = relFit, method = method,
                                           ...)
            if (verbose) {
                print(paste("Seed Network ", row, "/", n, sep = ""))
                if (!(verbose2 %in% "part")) {
                    print(toString(bitString))
                }
                print(paste(" - Score: ", fullScore, sep = ""))
                print("--------------------------------------------------")
                if (any(bitString != 0) & draw) {
                    plotDnf(model$reacID[which(bitString == 1)],
                            CNOlist = CNOlist)
                }
            }
            stop <- FALSE
            save.scores <- numeric()
            edges.changed <- numeric()
            edge.history <- character()
            counter <- 0
            while(!stop) {
                score <- computeScoreNemT1(CNOlist=CNOlist, model=model,
                                           bString=bitString, sizeFac=sizeFac,
                                           NAFac=NAFac, approach = approach,
                                           NEMlist = NEMlist,
                                           parameters=parameters, tellme = 0,
                                           relFit = relFit, method = method,
                                           ...)
                save.scores <-c(save.scores, score)
                scores <- numeric(bLength)
                sizes <- numeric(bLength)
                pValues <- numeric(bLength) + 1
                timeMark <- Sys.time()
                scoreThem <- function(i, CNOlist, model, simList = NULL,
                                      indexList, sizeFac, NAFac,
                                      approach = approach, NEMlist = NEMlist,
                                      parameters, tellme = 0, relFit = relFit,
                                      method = method, ...) {
                    if (!is.null(node)) {
                        if (length(grep(paste(node, collapse = "|"),
                                        model$reacID[i])) == 0) {
                            return(Inf)
                        } else {
                            bitStringTmp <- bitString
                            bitStringTmp[i] <- abs(1 - bitStringTmp[i])
                            redBstring <- reduceGraph(bitStringTmp, model,
                                                      CNOlist)
                            if (all(bitString == redBstring) & absorpII) {
                                bitStringTmp <- absorptionII(bitStringTmp,
                                                             model)
                            } else {
                                bitStringTmp <- redBstring
                            }
                            return(computeScoreNemT1(CNOlist=CNOlist,
                                                     model=model,
                                                     bString=bitStringTmp,
                                                     sizeFac=sizeFac,
                                                     NAFac=NAFac,
                                                     approach = approach,
                                                     NEMlist = NEMlist,
                                                     parameters=parameters,
                                                     tellme = 0,
                                                     relFit = relFit,
                                                     method = method,
                                                     ...))
                        }
                    } else {
                        bitStringTmp <- bitString
                        bitStringTmp[i] <- abs(1 - bitStringTmp[i])
                        redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
                        if (all(bitString == redBstring) & absorpII) {
                            bitStringTmp <- absorptionII(bitStringTmp, model)
                        } else {
                            bitStringTmp <- redBstring
                        }
                        if (sizeFac == 0) {
                            size <-
                                length(
                                    unlist(
                                        strsplit(
                                            model$reacID[
                                                      which(bitStringTmp == 1)],
                                            "\\+")))
                            return(c(computeScoreNemT1(CNOlist=CNOlist,
                                                       model=model,
                                                       bString=bitStringTmp,
                                                       sizeFac=sizeFac,
                                                       NAFac=NAFac,
                                                       approach = approach,
                                                       NEMlist = NEMlist,
                                                       parameters=parameters,
                                                       tellme = 0,
                                                       relFit = relFit,
                                                       method = method,
                                                       ...), size))
                        } else {
                            return(computeScoreNemT1(CNOlist=CNOlist,
                                                     model=model,
                                                     bString=bitStringTmp,
                                                     sizeFac=sizeFac,
                                                     NAFac=NAFac,
                                                     approach = approach,
                                                     NEMlist = NEMlist,
                                                     parameters=parameters,
                                                     tellme = 0,
                                                     relFit = relFit,
                                                     method = method,
                                                     ...))
                        }
                    }
                }
                edge.matrix <- as.matrix(seq_len(bLength))
                if (sum(prior != 0) > 0) {
                    edge.matrix <- as.matrix(edge.matrix[-which(prior != 0), ])
                }
                if (parallel2 == 1 & !is.null(parallel)) {
                    scores <- sfApply(edge.matrix, 1, scoreThem, CNOlist, model,
                                      simList = NULL, indexList, sizeFac, NAFac,
                                      approach = approach, NEMlist = NEMlist,
                                      parameters, tellme = 0, relFit = relFit,
                                      method = method)
                } else {
                    scores <- apply(edge.matrix, 1, scoreThem, CNOlist, model,
                                    simList = NULL, indexList, sizeFac, NAFac,
                                    approach = approach, NEMlist = NEMlist,
                                    parameters, tellme = 0, relFit = relFit,
                                    method = method)
                }
                if (sizeFac == 0) {
                    sizes <- scores[2, ]
                    scores <- scores[1, ]
                    size <-
                        length(
                            unlist(
                                strsplit(model$reacID[which(bitString == 1)],
                                         "\\+")))
                    check.size <- FALSE
                    if (any(scores == score) & all(scores >= score)) {
                        if (any(sizes[which(scores == score)] < size)) {
                            check.size <- TRUE
                        }
                    }
                } else {
                    check.size <- FALSE
                }
                if (sum(prior != 0) > 0) {
                    scores.tmp <- scores
                    scores <- numeric(length(model$reacID))
                    scores[which(prior == 0)] <- scores.tmp
                    scores[which(prior != 0)] <- Inf
                }
                timePassed <- as.numeric(Sys.time() - timeMark, unit = "secs")
                if (sum(scores < score) > 0 | check.size) {
                    if (check.size) { ##pValues[is.na(pValues)] <- 1
                        pValues <- scores # if you do not want to use pValues
                        topGates <- which(scores == score)
                        topScores <- scores[topGates]
                        topPvalues <- pValues[topGates]
                        minPvalue <- min(topPvalues)
                        topGate <- which(topPvalues == minPvalue)
                        topGate <- topGate[which(sizes[topGate] ==
                                                 min(sizes[topGate]))[1]]
                        topScore <- topScores[topGate]
                        whichGate <- which(sizes < size & scores == score)[1]
                        bitStringTmp <- bitString
                        bitStringTmp[whichGate] <-
                            abs(1 - bitStringTmp[whichGate])
                        redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
                        if (all(bitString == redBstring) & absorpII) {
                            bitStringTmp2 <- bitString
                            bitString <- absorptionII(bitStringTmp, model)
                            edges <- sum(bitString != bitStringTmp2)
                            edges.changed <- c(edges.changed, -edges)
                            if (verbose) {
                                deleted <-
                                    which(bitString == 0 & bitStringTmp == 1)
                                if (bitStringTmp[whichGate] == 0 &
                                    length(deleted) > 1) {
                                    deleted <-
                                        deleted[-which(deleted == whichGate)]
                                }
                                print(paste("Edges changed: ", edges, sep = ""))
                                print(paste("Deleted gates due to inverse ",
                                            "absorption: ",
                                            paste(model$reacID[deleted],
                                                  collapse = ", "), sep = ""))
                            }
                        } else {
                            bitStringTmp2 <- bitString
                            bitString <- redBstring
                            deleted <-
                                which(bitString == 0 & bitStringTmp2 == 1)
                            edges <- sum(bitString != bitStringTmp2)
                            edges.changed <- c(edges.changed, edges)
                            if (verbose & (length(deleted) > 1 |
                                           (length(deleted) > 0 &
                                            bitString[whichGate] == 1))) {
                                print(paste("Edges changed: ", edges, sep = ""))
                                print(paste("Deleted gates due to absorption: ",
                                            paste(model$reacID[deleted],
                                                  collapse = ", "), sep = ""))
                            }
                        }
                        bitStringsMem <- rbind(bitStringsMem, bitString)
                        if (bitString[whichGate] == 1) {
                            edge.history <- c(edge.history,
                                              model$reacID[whichGate])
                        }
                        if (verbose) {
                            counter <- counter + 1
                            print(paste("Iter step: ", counter, sep = ""))
                            if (bitString[whichGate] == 1) {
                                print(paste("Added gate ",
                                            model$reacID[whichGate], sep = ""))
                                if (!(verbose2 %in% "part")) {
                                    print(toString(bitString))
                                }
                                print(paste(" - Score: ", topScore, sep = ""))
                                print(paste(" - Iter_time: ", timePassed, " @ ",
                                            Sys.time(), sep = ""))
                                print("-------------------------
-------------------------")
                            } else {
                                print(paste("Deleted gate ",
                                            model$reacID[whichGate], sep = ""))
                                if (!(verbose2 %in% "part")) {
                                    print(toString(bitString))
                                }
                                print(paste(" - Score: ", topScore, sep = ""))
                                print(paste(" - Iter_time: ",
                                            timePassed, " @ ", Sys.time(),
                                            sep = ""))
                                print("---------------------------------------
-----------")
                            }
                            if (any(bitString != 0) & draw) {
                                plotDnf(model$reacID[which(bitString == 1)],
                                        CNOlist = CNOlist)
                            }
                        }
                        step.count <- step.count + 1
                        if (step.count >= max.steps) {
                            if (verbose) print("no more steps")
                            stop <- TRUE
                        }
                    } else {
                        ##pValues[is.na(pValues)] <- 1
                        pValues <- scores # if you do not want to use pValues
                        topGates <- which(scores < score)
                        topScores <- scores[topGates]
                        topPvalues <- pValues[topGates]
                        minPvalue <- min(topPvalues)
                        topGate <- which(topPvalues == minPvalue)[1]
                        topScore <- topScores[topGate]
                        whichGate <- topGates[topGate]
                        bitStringTmp <- bitString
                        bitStringTmp[whichGate] <-
                            abs(1 - bitStringTmp[whichGate])
                        redBstring <- reduceGraph(bitStringTmp, model, CNOlist)
                        if (all(bitString == redBstring) & absorpII) {
                            bitStringTmp2 <- bitString
                            bitString <- absorptionII(bitStringTmp, model)
                            edges <- sum(bitString != bitStringTmp2)
                            edges.changed <- c(edges.changed, -edges)
                            if (verbose) {
                                deleted <-
                                    which(bitString == 0 & bitStringTmp == 1)
                                if (bitStringTmp[whichGate] == 0 &
                                    length(deleted) > 1) {
                                    deleted <-
                                        deleted[-which(deleted == whichGate)]
                                }
                                print(paste("Edges changed: ", edges, sep = ""))
                                print(paste("Deleted gates due to inverse ",
                                            "absorption: ",
                                            paste(model$reacID[deleted],
                                                  collapse = ", "), sep = ""))
                            }
                        } else {
                            bitStringTmp2 <- bitString
                            bitString <- redBstring
                            deleted <-
                                which(bitString == 0 & bitStringTmp2 == 1)
                            edges <- sum(bitString != bitStringTmp2)
                            edges.changed <- c(edges.changed, edges)
                            if (verbose & (length(deleted) > 1 |
                                           (length(deleted) > 0 &
                                            bitString[whichGate] == 1))) {
                                print(paste("Edges changed: ", edges, sep = ""))
                                print(paste("Deleted gates due to absorption: ",
                                            paste(model$reacID[deleted],
                                                  collapse = ", "), sep = ""))
                            }
                        }
                        bitStringsMem <- rbind(bitStringsMem, bitString)
                        if (bitString[whichGate] == 1) {
                            edge.history <- c(edge.history,
                                              model$reacID[whichGate])
                        }
                        if (verbose) {
                            counter <- counter + 1
                            print(paste("Iter step: ", counter, sep = ""))
                            if (bitString[whichGate] == 1) {
                                print(paste("Added gate ",
                                            model$reacID[whichGate], sep = ""))
                                if (!(verbose2 %in% "part")) {
                                    print(toString(bitString))
                                }
                                print(paste(" - Score: ", topScore, sep = ""))
                                print(paste(" - Iter_time: ", timePassed, " @ ",
                                            Sys.time(), sep = ""))
                                print("---------------------------------------
-----------")
                            } else {
                                print(paste("Deleted gate ",
                                            model$reacID[whichGate], sep = ""))
                                if (!(verbose2 %in% "part")) {
                                    print(toString(bitString))
                                }
                                print(paste(" - Score: ", topScore, sep = ""))
                                print(paste(" - Iter_time: ", timePassed, " @ ",
                                            Sys.time(), sep = ""))
                                print("---------------------------------------
-----------")
                            }
                            if (any(bitString != 0) & draw) {
                                plotDnf(model$reacID[which(bitString == 1)],
                                        CNOlist = CNOlist)
                            }
                        }
                        step.count <- step.count + 1
                        if (step.count >= max.steps) {
                            if (verbose) print("no more steps")
                            stop <- TRUE
                        }
                    }
                } else {
                    if (verbose) print("no further improvement")
                    stop <- TRUE
                }
                if (Sys.time() - start.time > max.time) {
                    if (verbose) print("no more time")
                    stop <- TRUE
                }
            }
            return(list(A=score,B=bitString,row=row, C=save.scores,
                        D=edges.changed, E=edge.history))
        }
        if (!is.null(parallel) & parallel2 == 0) {
            bsTemp <- sfApply(as.matrix(seq_len(n)), 1, processSeed, bitStrings,
                              CNOlist, model, simList = NULL, indexList,
                              sizeFac, NAFac, approach = approach,
                              NEMlist = NEMlist, parameters, tellme = 0,
                              relFit = relFit, method = method,
                              max.steps = max.steps, node = node)
        } else {
            bsTemp <- apply(as.matrix(seq_len(n)), 1, processSeed, bitStrings,
                            CNOlist, model, simList = NULL, indexList, sizeFac,
                            NAFac, approach = approach, NEMlist = NEMlist,
                            parameters, tellme = 0, relFit = relFit,
                            method = method, max.steps = max.steps,
                            node = node)
        }
        save.scores <- list()
        edges.changed <- list()
        edge.history <- list()
        for (i in seq_len(length(bsTemp))) {
            bitStringsScores[bsTemp[[i]]$row] <- bsTemp[[i]]$A
            bitStrings[bsTemp[[i]]$row, ] <- bsTemp[[i]]$B
            save.scores[[i]] <- bsTemp[[i]]$C
            edges.changed[[i]] <- bsTemp[[i]]$D
            edge.history[[i]] <- bsTemp[[i]]$E
        }
        bitString <- bitStrings
        if (!is.null(parallel)) {
            sfStop()
        }
        bitString <- bitStrings # cbind(bitStrings, bitStringsScores)
        return(list(bStrings = bitString, scores = save.scores,
                    edges = edges.changed, history = edge.history))
    }
#' @noRd
makeDesignFull <-
    function(x, stimuli, inhibitors, batches = NULL, runs = NULL,
             method = "raw") {
        design0 <- makeDesign(x, stimuli, inhibitors, c(batches, runs))
        design <- design0

        ctrlsSum <- apply(design[, -grep(paste(c(runs, batches),
                                               collapse = "|"),
                                         colnames(design))], 1, sum)
        ctrlsSum <- which(ctrlsSum == 0)
        stimuliDesign <- design[, grep(paste(stimuli, collapse = "|"),
                                       colnames(design))]
        inhibitorsDesign <- design[, grep(paste(inhibitors, collapse = "|"),
                                          colnames(design))]
        if (length(stimuli) == 1) {
            stimuliDesign <- as.matrix(design[, grep(paste(stimuli,
                                                           collapse = "|"),
                                                     colnames(design))])
            colnames(stimuliDesign) <- stimuli
        }
        if (length(inhibitors) == 1) {
            inhibitorsDesign <- as.matrix(design[, grep(paste(inhibitors,
                                                              collapse = "|"),
                                                        colnames(design))])
            colnames(inhibitorsDesign) <- inhibitors
        }
        if (is.null(stimuli) == TRUE) {
            stimuliSum <- numeric(nrow(design))
        } else {
            if (length(stimuli) == 1) {
                stimuliSum <- stimuliDesign
            } else {
                stimuliSum <- apply(stimuliDesign, 1, sum)
            }
        }
        if (is.null(inhibitors) == TRUE) {
            inhibitorsSum <- numeric(nrow(design))
        } else {
            if (length(inhibitors) == 1) {
                inhibitorsSum <- inhibitorsDesign
            } else {
                inhibitorsSum <- apply(inhibitorsDesign, 1, sum)
            }
        }
        cuesSum <- apply(design[, grep(paste(c(stimuli, inhibitors),
                                             collapse = "|"),
                                       colnames(design))], 1, sum)
        maxStim <- max(stimuliSum)
        maxKd <- max(inhibitorsSum)
        maxCue <- max(cuesSum)

        grepCtrl <- which(cuesSum == 0)
        grepStims <- intersect(which(stimuliSum != 0),
                               which(inhibitorsSum == 0))
        grepKds <- intersect(which(stimuliSum == 0),
                             which(inhibitorsSum != 0))
        grepStimsKds <- intersect(which(stimuliSum != 0),
                                  which(inhibitorsSum != 0))

        design <- numeric()
        designNames <- character()
        design <- rep(0, ncol(x))
        design[grepCtrl] <- 1
        designNames <- "Ctrl"

        for (i in grepStims) {
            stimNames <- paste(sort(names(which(stimuliDesign[i, ] >= 1))),
                               collapse = "_")
            if (stimNames %in% designNames) {
                design[i, which(designNames %in% stimNames)] <- 1
            } else {
                design <- cbind(design, rep(0, ncol(x)))
                designNames <- c(designNames, stimNames)
                design[i, which(designNames %in% stimNames)] <- 1
            }
        }

        for (i in grepKds) {
            stimNames <- paste(sort(names(which(inhibitorsDesign[i, ] >= 1))),
                               collapse = "_")
            if (stimNames %in% designNames) {
                design[i, which(designNames %in% stimNames)] <- 1
            } else {
                design <- cbind(design, rep(0, ncol(x)))
                designNames <- c(designNames, stimNames)
                design[i, which(designNames %in% stimNames)] <- 1
            }
        }

        for (i in grepStimsKds) {
            stimNames <- paste(c(sort(names(which(inhibitorsDesign[i, ] >= 1))),
                                 sort(names(which(stimuliDesign[i, ] >= 1)))),
                               collapse = "_")
            if (stimNames %in% designNames) {
                design[i, which(designNames %in% stimNames)] <- 1
            } else {
                design <- cbind(design, rep(0, ncol(x)))
                designNames <- c(designNames, stimNames)
                design[i, which(designNames %in% stimNames)] <- 1
            }
        }

        if (!is.null(batches))  {
            for (i in batches) {
                if (!is.null(runs)) {
                    for (j in runs) {
                        tmp <- numeric(ncol(x))
                        tmp[intersect(grep(i, colnames(x)),
                                      grep(j, colnames(x)))] <- 1
                        if (sum(tmp) != 0) {
                            design <- cbind(design, tmp)
                            designNames <- c(designNames, paste(sort(c(i, j)),
                                                                collapse = "_"))
                        }
                    }
                }
            }
        }
        colnames(design) <- designNames
        if (method %in% "inter") {
            for (i in c(stimuli, inhibitors)) {
                design[, i] <- design0[, i]
            }
        }
        return(design)
    }
#' @noRd
makeDesign <-
    function(x, stimuli, inhibitors, batches = NULL, runs = NULL) {
        design <- numeric()
        designNames <- character()
        for (i in c(stimuli, inhibitors, batches, runs)) {
            tmp <- numeric(ncol(x))
            tmp[grep(i, colnames(x))] <- 1
            if (sum(tmp) != 0) {
                design <- cbind(design, tmp)
                designNames <- c(designNames, i)
            }
        }
        colnames(design) <- designNames
        return(design)
    }
#' @noRd
#' @import snowfall
myGsea <-
    function(testList, goList, parallel = NULL, adjust.method = "FDR",
             conservative = TRUE) {
        if (!is.null(parallel)) {
            if (is.list(parallel)) {
                if (length(parallel[[1]]) != length(parallel[[2]])) {
                    stop("The nodes (second list object in parallel) and the
number of cores used on every node (first list object in parallel) must be
the same.") }
                hosts <- character()
                for (i in seq_len(length(parallel[[1]]))) {
                    hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
                }
                hosts <- as.list(hosts)
                sfInit(parallel=TRUE, socketHosts=hosts)
            } else {
                sfInit(parallel=TRUE, cpus=parallel)
            }
            sfExport("testList", "goList")
        }

        startTime <- Sys.time()

        if (conservative) {
            Complete <- unique(intersect(unlist(testList),unlist(goList)))
        } else {
            Complete <- unique(c(unlist(testList),unlist(goList)))
        }
        N <- length(Complete)

        checkGeneCluster <- function(j,i) {

            resDf <- data.frame()
            genelist <- unique(intersect(testList[[i]],Complete))

            if (!is.null(genelist)) {
                targetlist <- intersect(goList[[j]],Complete)
                n <- length(targetlist)
                D <- length(genelist)
                AinB <- length(intersect(targetlist, genelist))
                fisherMat <- cbind(c(AinB,D-AinB),c(n-AinB,N-n-D+AinB))
                fisherTest <- fisher.test(fisherMat, alternative = "greater")
                ABpval <- fisherTest$p.value
                resDf <- rbind(resDf, data.frame(test = names(testList)[i],
                                                 go = names(goList)[j],
                                                 pval = ABpval, overlap = AinB,
                                                 targetset = n, testset = D))
            }

            return(resDf)
        }

        dfList <- NULL
        if (is.null(parallel)) {
            for (i in seq_len(length(testList))) {
                cat(".")
                dfList <- c(dfList, lapply(seq_len(length(goList)),
                                           checkGeneCluster, i))
            }
        } else {
            for (i in seq_len(length(testList))) {
                dfList <- c(dfList, sfLapply(seq_len(length(goList)),
                                             checkGeneCluster, i))
            }
        }

        cat("\n")
        print(Sys.time() - startTime)
        print("all hypergeometric tests (=fisher tests with 'greater')
 completed...")
        print("preparing data...")

        bigDf <- do.call("rbind", dfList)

        bigDf <- bigDf[order(bigDf[, 3]), ]

        if (adjust.method %in% "FDR") {

            qvals <- (nrow(bigDf)*bigDf[, 3])/(seq_len(nrow(bigDf)))

            for (i in (length(qvals)-1):1) {

                qvals[i] <- min(qvals[i], qvals[i+1])

            }

            qvals[which(qvals > 1)] <- 1

        } else {

            qvals <- p.adjust(bigDf[, 3], method = adjust.method)

        }

        bigDf <- cbind(bigDf, qvals = qvals)

        bigDf <- bigDf[order(bigDf[, 7]), ]

        if (!is.null(parallel)) {
            sfStop()
        }

        return(bigDf)
    }
#' @noRd
#' @import cluster
pamNorm <-
    function(x, method = "euclidean") {
        if (is.matrix(x)) {
            for (i in seq_len(nrow(x))) {
                clust <- pam(x[i, ], 2)
                clust.nums <- clust$clustering
                dist <- dist(x[i, ], method = method)
                sil <- silhouette(clust.nums, dist=dist)
                max <- which(x[i, ] == max(x[i, ]))
                cluster.one <- sil[max, 1]
                cluster.zero <- c(1,2)[-cluster.one]
                min <- which(x[i, ] == min(x[i, ]))
                sil[which(sil[, 1] == cluster.zero), 3] <-
                    sil[which(sil[, 1] == cluster.zero), 3]*(-1)
                sil[, 3] <- (sil[, 3] + 1)/2
                sil[which(x[i, ] >= x[i, which(sil[, 3] ==
                                               max(sil[, 3]))[1]]), 3] <-
                    sil[which(x[i, ] == x[i, which(sil[, 3] ==
                                                   max(sil[, 3]))[1]])[1], 3]
                sil[which(x[i, ] <= x[i, which(sil[, 3] ==
                                               min(sil[, 3]))[1]]), 3] <-
                    sil[which(x[i, ] == x[i, which(sil[, 3] ==
                                                   min(sil[, 3]))[1]])[1], 3]
                sil[, 3] <- sil[, 3] - min(sil[, 3])
                sil[, 3] <- sil[, 3]/max(sil[, 3])
                x[i, ] <- sil[, 3]
            }
        } else {
            clust <- pam(x, 2)
            clust.nums <- clust$clustering
            dist <- dist(x, method = method)
            sil <- silhouette(clust.nums, dist=dist)
            max <- which(x == max(x))
            cluster.one <- sil[max, 1]
            cluster.zero <- c(1,2)[-cluster.one]
            min <- which(x == min(x))
            sil[which(sil[, 1] == cluster.zero), 3] <-
                sil[which(sil[, 1] == cluster.zero), 3]*(-1)
            sil[, 3] <- (sil[, 3] + 1)/2
            sil[which(x >= x[which(sil[, 3] == max(sil[, 3]))[1]]), 3] <-
                sil[which(x == x[which(sil[, 3] == max(sil[, 3]))[1]])[1], 3]
            sil[which(x <= x[which(sil[, 3] == min(sil[, 3]))[1]]), 3] <-
                sil[which(x == x[which(sil[, 3] == min(sil[, 3]))[1]])[1], 3]
            sil[, 3] <- sil[, 3] - min(sil[, 3])
            sil[, 3] <- sil[, 3]/max(sil[, 3])
            x <- sil[, 3]
        }
        return(x)
    }
#' @noRd
removeCycles <-
    function(bString, model, dnf = NULL) {
        if (is.null(dnf)) {
            if (any(bString != 0)) {
                graph <- model$reacID[which(bString == 1)]
                adjmat <- abs(dnf2adj(graph))
                get.order <- apply(adjmat, 2, sum)
                adjmat <- adjmat[order(get.order), order(get.order)]
                cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1,
                                arr.ind = TRUE)
                while(any(adjmat[lower.tri(adjmat)] == 1) &
                      dim(cycles)[1] > 0) {
                    for (i in seq_len(nrow(cycles))) {
                        bad.cycles <-
                            grep(paste(".*", rownames(adjmat)[cycles[i, 1]],
                                       ".*=", colnames(adjmat)[cycles[i, 2]],
                                       sep = ""), model$reacID)
                        if (length(bad.cycles) > 0 &
                            any(bString[bad.cycles] == 1)) {
                            bString[bad.cycles] <- 0
                            graph <- model$reacID[which(bString == 1)]
                            break()
                        }
                    }
                    adjmat <- abs(dnf2adj(graph))
                    get.order <- apply(adjmat, 2, sum)
                    adjmat <- adjmat[order(get.order), order(get.order)]
                    cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1,
                                    arr.ind = TRUE)
                }
                return(bString)
            } else {
                return(bString)
            }
        } else {
            graph <- dnf
            adjmat <- abs(dnf2adj(graph))
            get.order <- apply(adjmat, 2, sum)
            adjmat <- adjmat[order(get.order), order(get.order)]
            cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1,
                            arr.ind = TRUE)
            while(any(adjmat[lower.tri(adjmat)] == 1) & dim(cycles)[1] > 0) {
                for (i in seq_len(nrow(cycles))) {
                    bad.cycles <-
                        grep(paste(".*",
                                   rownames(adjmat)[cycles[i, 1]],
                                   ".*=",
                                   colnames(adjmat)[cycles[i, 2]], sep = ""),
                             graph)
                    if (length(bad.cycles) > 0) {
                        graph <- graph[-bad.cycles]
                        break()
                    }
                }
                adjmat <- abs(dnf2adj(graph))
                get.order <- apply(adjmat, 2, sum)
                adjmat <- adjmat[order(get.order), order(get.order)]
                cycles <- which(lower.tri(adjmat) == TRUE & adjmat == 1,
                                arr.ind = TRUE)
            }
            return(graph)
        }
    }
#' @noRd
simpleNorm <-
    function(x) {
        if (is.matrix(x) == FALSE) {
            genes.min <- min(x)
            x <- (x - genes.min)
            genes.max <- max(x)
            x <- x/genes.max
        } else {
            genes.min <- apply(x, 1, min)
            x <- (x - genes.min)
            genes.max <- apply(x, 1, max)
            x <- x/genes.max
        }
        x[is.na(x)] <- 0
        return(x)
    }
#' @noRd
#' @import matrixStats
simulateDnf <-
    function(dnf, stimuli = NULL, inhibitors = NULL) {
        getStateDnf <- function(node, signalStates, graph, children = NULL) {
            graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
            if (length(graphCut) == 0) {
                signalStates[, node] <- 0
            } else {
                sop <- numeric(nrow(signalStates))
                children2 <- gsub("!", "", children)
                for (i in graphCut) {
                    parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
                    pob <- rep(1, nrow(signalStates))
                    for (j in parents) {
                        j2 <- gsub("!", "", j)
                        if (sum(is.na(signalStates[, j2]) == TRUE) ==
                            length(signalStates[, j2])) {
                            if (j %in% j2) {
                                node2 <- node
                                add1 <- 0
                            } else {
                                node2 <- paste("!", node, sep = "")
                                add1 <- 1
                            }
                            if (j2 %in% children2) {
                                subGraph <-
                                    graph[-grep(paste(".*=", node, "|.*", j2,
                                                      ".*=.*", sep = ""),
                                                graph)]
                                signalStatesTmp <-
                                    getStateDnf(node = j2,
                                                signalStates = signalStates,
                                                graph = subGraph,
                                                children = NULL)
                                if (
                                (length(
                                    grep(
                                        "!",
                                        children[which(children2 %in%
                                                       j2):
                                                 length(
                                                     children2)]))+add1)/2 !=
                                ceiling(
                                (
                                    length(
                                        grep("!",
                                             children[which(children2 %in%
                                                            j2):
                                                      length(children2)]))+
                                    add1)/2)) {
                                    ## negative feedback loop calculation does
                                    ## not seem to be general enough and also
                                    ## not feasible:
                                } else {
                                }
                                if (add1 == 0) {
                                    pobMult <- signalStatesTmp[, j2]
                                } else {
                                    pobMult <- add1 - signalStatesTmp[, j2]
                                }
                            } else {
                                signalStates <-
                                    getStateDnf(node = j2,
                                                signalStates = signalStates,
                                                graph = graph,
                                                children = unique(c(children,
                                                                    node2)))
                                if (add1 == 0) {
                                    pobMult <- signalStates[, j2]
                                } else {
                                    pobMult <- add1 - signalStates[, j2]
                                }
                            }
                            pob <- pob*pobMult
                        } else {
                            if (j %in% j2) {
                                add1 <- 0
                            } else {
                                add1 <- 1
                            }
                            if (add1 == 0) {
                                pobMult <- signalStates[, j2]
                            } else {
                                pobMult <- add1 - signalStates[, j2]
                            }
                            pob <- pob*pobMult
                        }
                        if (max(pob, na.rm = TRUE) == 0) { break() }
                    }
                    sop <- sop + pob
                    if (min(sop, na.rm = TRUE) > 0) { break() }
                }
                sop[sop > 0] <- 1
                if (node %in% inhibitors) {
                    sop <- sop*0
                }
                if (node %in% stimuli) {
                    sop <- max(sop, 1)
                }
                signalStates[, node] <- sop
            }
            return(signalStates)
        }
        signals <-
            unique(gsub("!", "", unlist(strsplit(unlist(strsplit(dnf, "=")),
                                                 "\\+"))))
        graph <- dnf
        signalStates <- matrix(NA, nrow = 1, ncol = length(signals))
        rownames(signalStates) <- paste(c("stimuli:", stimuli, "inhibitors:",
                                          inhibitors), collapse = " ")
        colnames(signalStates) <- signals
        signalStates[which(signals %in% stimuli)] <- 1
        for (k in signals) {
            if (is.na(signalStates[, k]) == TRUE) {
                signalStates <- getStateDnf(node = k,
                                            signalStates = signalStates,
                                            graph = graph, children = NULL)
            }
        }
        namestmp <- colnames(signalStates)
        signalStates <- as.vector(signalStates)
        names(signalStates) <- namestmp
        return(signalStates = signalStates)
    }
#' @noRd
#' @import
#' matrixStats
#' @importFrom Biobase rowMin rowMax
simulateStatesRecursiveAdd <-
    function(CNOlist, model, bString, NEMlist = NULL) {
        getStateAdd <- function(CNOlist, node, signalStates, graph,
                                children = NULL, NEMlist = NULL) {
            graphCut <- graph[grep(paste("=", node, sep = ""), graph)]
            if (length(graphCut) == 0) {
                if (node %in% colnames(CNOlist@inhibitors)) {
                    signalStates[, node] <- 0 - CNOlist@inhibitors[, node]
                } else {
                    signalStates[, node] <- 0
                }
            } else {
                sop <- numeric(nrow(signalStates)) - 1
                children2 <- gsub("!", "", children)
                for (i in graphCut) {
                    parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
                    pob <- rep(1, nrow(signalStates))
                    for (j in parents) {
                        j2 <- gsub("!", "", j)
                        if (sum(is.na(signalStates[, j2]) == TRUE) ==
                            length(signalStates[, j2])) {
                            if (j2 %in% children2) {
                                if (j2 %in% j) {
                                    node2 <- node
                                    add1 <- 0
                                } else {
                                    node2 <- paste("!", node, sep = "")
                                    add1 <- 1
                                }
                                if (
                                (length(
                                    grep(
                                        "!",
                                        children[
                                            which(
                                                children2 %in%
                                                j2):length(children2)]))+
                                 add1)/2 !=
                                ceiling((
                                    length(
                                        grep(
                                            "!",
                                            children[
                                                which(
                                                    children2 %in% j2):
                                                length(
                                                    children2)]))+add1)/2)) {
                                    subGraph <-
                                        graph[
                                            -grep(paste(".*",
                                                        node,
                                                        ".*=.*",
                                                        children2[
                                                            length(children2)],
                                                        sep = ""), graph)]
                                    subResult <-
                                        getStateAdd(CNOlist = CNOlist,
                                                    node = node,
                                                    signalStates = signalStates,
                                                    graph = subGraph,
                                                    children = NULL, NEMlist)
                                    pobMult <- subResult[, node]
                                    subResult <- signalStates
                                    subResult[, node] <- pobMult
                                    subGraph2 <- graph[-which(graph %in% i)]
                                    subResult2 <-
                                        getStateAdd(CNOlist = CNOlist,
                                                    node = j2,
                                                    signalStates = subResult,
                                                    graph = subGraph2,
                                                    children = NULL, NEMlist)
                                    if (add1 == 0) {
                                        pobMult2 <- subResult2[, j2]
                                    } else {
                                        pobMult2 <- add1 - subResult2[, j2]
                                    }
                                    pobMult2[pobMult2 == 2] <- 1

                                    pobNA <- numeric(length(pob))
                                    pobNA[is.na(pob)] <- 1
                                    pobNA[is.na(pobMult)] <- 1
                                    pobNA[which(pobMult != pobMult2)] <- 1
                                    pobMult[is.na(pobMult)] <- 1
                                    pobMult[which(pobMult != pobMult2)] <- 1
                                    pob[is.na(pob)] <- 1

                                    ##pobMult[pobMult == -1] <- 0
                                    pob <- rowMin(cbind(pob,pobMult))

                                    pobNA[which(pob == 0)] <- 0
                                    pob[which(pobNA > 0)] <- NA
                                } else {
                                    signalStatesTemp <- signalStates
                                    if (!(j %in% j2)) {
                                        if (j2 %in%
                                            colnames(CNOlist@inhibitors)) {
                                            signalStatesTemp[, node] <-
                                                1 - CNOlist@inhibitors[, j2]
                                        } else {
                                            signalStatesTemp[, node] <- 1
                                        }
                                    } else {
                                        if (j2 %in%
                                            colnames(CNOlist@inhibitors)) {
                                            signalStatesTemp[, node] <-
                                                0 - CNOlist@inhibitors[, j2]
                                        } else {
                                            signalStatesTemp[, node] <- 0
                                        }
                                    }
                                    subGraph <- graph[-which(graph %in% i)]
                                    subResult <-
                                        getStateAdd(CNOlist = CNOlist,
                                                    node = j2,
                                                    signalStates =
                                                        signalStatesTemp,
                                                    graph = subGraph,
                                                    children = NULL, NEMlist)

                                    if (add1 == 0) {
                                        pobMult <- subResult[, j2]
                                    } else {
                                        pobMult <- add1 - subResult[, j2]
                                    }
                                    pobMult[pobMult == 2] <- 1

                                    pobNA <- numeric(length(pob))
                                    pobNA[is.na(pob)] <- 1
                                    pobNA[is.na(pobMult)] <- 1
                                    pobMult[is.na(pobMult)] <- 1
                                    pob[is.na(pob)] <- 1

                                    ##pobMult[pobMult == -1] <- 0
                                    pob <- rowMin(cbind(pob,pobMult))

                                    pobNA[which(pob == 0)] <- 0
                                    pob[which(pobNA > 0)] <- NA
                                }
                            } else {
                                if (j %in% j2) {
                                    node2 <- node
                                    add1 <- 0
                                } else {
                                    node2 <- paste("!", node, sep = "")
                                    add1 <- 1
                                }
                                signalStates <-
                                    getStateAdd(CNOlist = CNOlist, node = j2,
                                                signalStates = signalStates,
                                                graph = graph,
                                                children =
                                                    unique(c(children, node2)),
                                                NEMlist)

                                if (add1 == 0) {
                                    pobMult <- signalStates[, j2]
                                } else {
                                    pobMult <- add1 - signalStates[, j2]
                                }
                                pobMult[pobMult == 2] <- 1

                                pobNA <- numeric(length(pob))
                                pobNA[is.na(pob)] <- 1
                                pobNA[is.na(pobMult)] <- 1
                                pobMult[is.na(pobMult)] <- 1
                                pob[is.na(pob)] <- 1

                                ##pobMult[pobMult == -1] <- 0
                                pob <- rowMin(cbind(pob,pobMult))

                                pobNA[which(pob == 0)] <- 0
                                pob[which(pobNA > 0)] <- NA
                            }
                        } else {
                            if (j %in% j2) {
                                add1 <- 0
                            } else {
                                add1 <- 1
                            }

                            if (add1 == 0) {
                                pobMult <- signalStates[, j2]
                            } else {
                                pobMult <- add1 - signalStates[, j2]
                            }
                            pobMult[pobMult == 2] <- 1

                            pobNA <- numeric(length(pob))
                            pobNA[is.na(pob)] <- 1
                            pobNA[is.na(pobMult)] <- 1
                            pobMult[is.na(pobMult)] <- 1
                            pob[is.na(pob)] <- 1
                            pob <- rowMin(cbind(pob,pobMult))

                            pobNA[which(pob == 0)] <- 0
                            pob[which(pobNA > 0)] <- NA
                        }
                        if (max(pob, na.rm = TRUE) == 0) { break() }
                    }
                    pobNA <- numeric(length(pob))
                    pobNA[is.na(pob)] <- 1
                    pobNA[is.na(sop)] <- 1
                    pob[is.na(pob)] <- 0
                    sop[is.na(sop)] <- 0
                    sop <- rowMax(cbind(sop,pob))
                    pobNA[which(sop > 0)] <- 0
                    sop[which(pobNA > 0)] <- NA
                    if (min(sop, na.rm = TRUE) > 0) { break() }
                }
                ##sop[sop > 0] <- 1
                if (node %in% colnames(CNOlist@inhibitors)) {
                    sop <- sop - CNOlist@inhibitors[, node]
                }
                if (node %in% colnames(CNOlist@stimuli)) {
                    sop <- sop + CNOlist@stimuli[, node]
                }
                signalStates[, node] <- sop
            }
            return(signalStates)
        }
        bString <- reduceGraph(bString, CNOlist, model)
        stimuli <- colnames(CNOlist@stimuli)
        inhibitors <-
            c(colnames(CNOlist@inhibitors),
              model$namesSpecies[-which(model$namesSpecies %in%
                                        c(stimuli,
                                          colnames(CNOlist@inhibitors)))])
        graph <- model$reacID[which(bString == 1)]
        stimuliStates <- CNOlist@stimuli
        if (!is.null(NEMlist$signalStates)) {
            signalStates <- NEMlist$signalStates
        } else {
            signalStates <- matrix(NA, nrow = nrow(CNOlist@signals[[2]]),
                                   ncol = length(inhibitors))
            rownames(signalStates) <- rownames(CNOlist@signals[[2]])
            colnames(signalStates) <- inhibitors
            signalStates <- cbind(stimuliStates, signalStates)
        }
        for (k in inhibitors) {
            if (sum(is.na(signalStates[, k]) == TRUE) ==
                length(signalStates[, k])) {
                signalStates <- getStateAdd(CNOlist = CNOlist, node = k,
                                            signalStates = signalStates,
                                            graph = graph, children = NULL,
                                            NEMlist)
            }
        }
        return(signalStates)
    }
#' @noRd
smoothMatrix <-
    function(M, n=1, direction = 0, torus = FALSE) {
        Msmooth <- M
        if (n > 0) {
            for (i in seq_len(n)) {

                cat('\r', i)
                flush.console()

                Mtmp <- Msmooth
                M1 <- M2 <- M3 <- M4 <- M5 <- M6 <- M7 <- M8 <- M*0
                if (torus) {
                    if (any(direction %in% c(0,1))) {
                        M1 <- cbind(Msmooth[, ncol(M)],
                                    Msmooth[, seq_len((ncol(M)-1))])
                        M5 <- cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])
                    }
                    if (any(direction %in% c(0,2))) {
                        M3 <- rbind(Msmooth[nrow(M), ],
                                    Msmooth[seq_len((nrow(M)-1)), ])
                        M7 <- rbind(Msmooth[2:nrow(M), ], Msmooth[1, ])
                    }
                    if (any(direction %in% c(0,3))) {
                        M2 <- rbind(cbind(Msmooth[, ncol(M)],
                                          Msmooth[
                                            , seq_len((ncol(M)-1))])[nrow(M), ],
                                    cbind(Msmooth[, ncol(M)],
                                          Msmooth[
                                            , seq_len((ncol(M)-1))])[
                                        seq_len((nrow(M)-1)), ])
                        M6 <- rbind(cbind(Msmooth[, 2:ncol(M)],
                                          Msmooth[, 1])[2:nrow(M), ],
                                    cbind(Msmooth[, 2:ncol(M)],
                                          Msmooth[, 1])[1, ])
                    }
                    if (any(direction %in% c(0,4))) {
                        M4 <- rbind(cbind(Msmooth[, 2:ncol(M)],
                                          Msmooth[, 1])[nrow(M), ],
                                    cbind(Msmooth[, 2:ncol(M)],
                                          Msmooth[, 1])[seq_len((nrow(M)-1)), ])
                        M8 <- cbind(rbind(Msmooth[2:nrow(M), ],
                                          Msmooth[1, ])[, ncol(M)],
                                    rbind(Msmooth[2:nrow(M), ],
                                          Msmooth[1, ])[, seq_len((ncol(M)-1))])
                    }
                } else {
                    if (any(direction %in% c(0,1))) {
                        M1 <- cbind(0, Msmooth[, seq_len((ncol(M)-1))])
                        M5 <- cbind(Msmooth[, 2:ncol(M)], 0)
                    }
                    if (any(direction %in% c(0,2))) {
                        M3 <- rbind(0, Msmooth[seq_len((nrow(M)-1)), ])
                        M7 <- rbind(Msmooth[2:nrow(M), ], 0)
                    }
                    if (any(direction %in% c(0,3))) {
                        M2 <-
                            rbind(0,
                                  cbind(0,
                                        Msmooth[
                                          , seq_len((ncol(M)-1))])[
                                      seq_len((nrow(M)-1)), ])
                        M6 <- rbind(cbind(Msmooth[, 2:ncol(M)], 0)[2:nrow(M), ]
                                  , 0)
                    }
                    if (any(direction %in% c(0,4))) {
                        M4 <- rbind(0, cbind(Msmooth[, 2:ncol(M)], 0)[
                                           seq_len((nrow(M)-1)), ])
                        M8 <- cbind(0, rbind(Msmooth[2:nrow(M), ], 0)[
                                         , seq_len((ncol(M)-1))])
                    }
                }
                denom <- matrix(9, nrow(M), ncol(M))
                Msmooth <- Mtmp+M1+M2+M3+M4+M5+M6+M7+M8
                Msmooth <- Msmooth/denom
                if (all(Mtmp == Msmooth)) {
                    print("convergence")
                    break()
                }
            }
        }
        return(Msmooth)
    }
#' @noRd
createCube <- function(n=3, m=n) {
    if (m > n) { m <- n }
    n2 <- FALSE
    m2 <- FALSE
    if (!(round(n/2) == n/2)) {
        n2 <- TRUE
        n <- n + 1
    }
    if (!(round(m/2) == m/2)) {
        m2 <- TRUE
        m <- m + 1
    }
    e <- matrix(0, m, n)
    e[1, ] <- sample(c(0,1), n, replace = TRUE)
    for (i in (seq_len((m/2))*2)) {
        e[i, ] <- 1 - e[i-1, ]
        if (i >= m) { break() }
        for (j in seq_len(i)) {
            e[i+1, ((j-1)*(n/i)+1):(j*(n/i))] <- e[j, ((j-1)*(n/i)+1):(j*(n/i))]
        }
    }
    if (n2) {
        e <- e[, -1]
    }
    if (m2) {
        e <- e[-nrow(e), ]
    }
    return(e)
}
