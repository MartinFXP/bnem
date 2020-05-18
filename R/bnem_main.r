#' BCR perturbation reproduction
#'
#' Produce the application data from the BCR paper
#' of Pirkl, et al., 2016, Bioinformatics. Raw data is available at
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68761
#' @param path path to the CEL.gz files
#' @param combsign if TRUE includes all covariates in ComBat analysis
#' @author Martin Pirkl
#' @return list with the full foldchanges and epxression matrix,
#' a reduced foldchange matrix and the design matrix for the computations
#' @export
#' @importFrom vsn justvsn
#' @importFrom limma lmFit eBayes makeContrasts contrasts.fit
#' @importFrom sva ComBat
#' @importFrom affy rma ReadAffy list.celfiles
#' @importFrom Biobase exprs
#' @examples
#'\dontrun{
#' processDataBCR()
#' }
#' data(bcr)
processDataBCR <- function(path = "", combsign = FALSE) {
    wd <- getwd()
    setwd(path)
    fns <- affy::list.celfiles()
    ab <- affy::ReadAffy(filenames=fns)
    setwd(wd)
    ev <- vsn::vsnrma(ab)
    data <- exprs(ev)
    colnames(data) <- gsub("GSM16808[0-9][0-9]_|\\.CEL\\.gz", "",
                          colnames(data))

    batch <- gsub(".*_", "", colnames(data))
    batch[which(batch == "Batch")] <- "Batch3"
    colnames(data) <- gsub("SP600125", "Jnk", colnames(data))
    colnames(data) <- gsub("SB203580", "p38", colnames(data))
    colnames(data) <- gsub("10058F4", "Myc", colnames(data))
    colnames(data) <- gsub("Ko", "KO", colnames(data))
    colnames(data) <- gsub("Ly294002", "LY294002", colnames(data))
    colnames(data) <- gsub("LY294002", "Pi3k", colnames(data))
    colnames(data) <- gsub("U0126", "Erk", colnames(data))
    colnames(data) <- gsub("DMSO", "Ctrl", colnames(data))
    colnames(data) <- gsub("IKK2", "Ikk2", colnames(data))
    colnames(data) <- gsub("TAK1", "Tak1", colnames(data))
    vars <- unique(unlist(strsplit(gsub("_Batch.*", "", colnames(data)),
                                   "_")))
    vars <- sort(vars[-grep("\\+", vars)])
    vars <- vars[-which(vars %in% c("KO"))]
    design <- matrix(0, ncol(data), length(vars))
    colnames(design) <- vars
    rownames(design) <- colnames(data)
    for (i in vars) {
        design[grep(i, colnames(data)), i] <- 1
    }
    combos <- NULL
    for (i in which(apply(design, 1, sum) > 1)) {
        combos <- c(combos,
                    paste(sort(colnames(design)[which(design[i, ] > 0)]),
                          collapse = "_"))
    }
    combos <- sort(unique(combos))
    design2 <- design
    design2[which(apply(design, 1, sum) > 1), ] <- 0
    for (i in combos) {
        comb <- unlist(strsplit(i, "_"))
        tmp2 <- numeric(nrow(design))
        tmp2[which(apply(design[, comb], 1, sum) == length(comb) &
                  apply(design, 1, sum) == length(comb))] <- 1
        design2 <- cbind(design2, tmp2)
        colnames(design2)[ncol(design2)] <- i
    }
    design2 <- design2[, -which(colnames(design2) %in% "BCR")]
    colnames(design2)[which(colnames(design2) %in% "BCR_Ctrl")] <- "BCR"

    if (combsign) {
        dataCB <- sva::ComBat(data, batch,
                              design2[, -grep("Ctrl", colnames(design2))])
    } else {
        dataCB <- sva::ComBat(data, batch)
    }

    fit <- limma::lmFit(dataCB, design2)

    fc <- matrix(0, nrow(dataCB), (ncol(design2)-2)*2 + 1)
    colnames(fc) <- seq_len(ncol(fc))
    contmat <- limma::makeContrasts(Ctrl_vs_BCR="BCR-Ctrl", levels=design2)
    fit2 <- limma::contrasts.fit(fit, contmat)
    fit2 <- limma::eBayes(fit2)
    fc[, 1] <- fit2$coefficients
    colnames(fc)[1] <- "Ctrl_vs_BCR"
    colnames(fc)[-1] <- c(paste("Ctrl_vs",
                                colnames(design2)[-which(colnames(design2) %in%
                                                         c("Ctrl", "BCR"))],
                                sep = "_"),
                          paste("BCR_vs",
                                colnames(design2)[-which(colnames(design2) %in%
                                                         c("Ctrl", "BCR"))],
                                sep = "_"))
    for (i in colnames(design2)) {
        if (i %in% c("Ctrl", "BCR")) { next() }
        contmat <- contmat*0
        contmat[which(colnames(design2) %in% "Ctrl")] <- -1
        contmat[i, ] <- 1
        fit2 <- limma::contrasts.fit(fit, contmat)
        fit2 <- limma::eBayes(fit2)
        fc[, paste0("Ctrl_vs_", i)] <- fit2$coefficients
        contmat <- contmat*0
        contmat[which(colnames(design2) %in% "BCR")] <- -1
        contmat[i, ] <- 1
        fit2 <- limma::contrasts.fit(fit, contmat)
        fit2 <- limma::eBayes(fit2)
        fc[, paste0("BCR_vs_", i)] <- fit2$coefficients
    }
    targets <- paste("BCR_vs_BCR",
                                colnames(design2)[-which(colnames(design2) %in%
                                                         c("DMSO", "BCR"))],
                                sep = "_")
    targets <- targets[-grep("Myc|LY294|U0126|Vivit|BCR_BCR|BCR_Ctrl", targets)]
    fc2 <- fc[, c("Ctrl_vs_BCR", targets)]
    rownames(fc) <- rownames(data)

    fc2 <- fc2[which(abs(fc2[, "Ctrl_vs_BCR"]) > 1 &
                   apply(abs(fc2[, -which(colnames(fc2) %in%
                                         "Ctrl_vs_BCR")]), 1, max) >
                   log2(1.5)), ]
    fci <- fc2[, -which(colnames(fc2) %in%
                        "Ctrl_vs_BCR")]*sign(fc2[, "Ctrl_vs_BCR"])
    argl <- apply(fci, 1, min)
    fc2 <- fc2[which(argl < 0), ]

    bcr <- list(exprs = fit$coefficients%*%t(design2), fc = fc2, full = fc,
                design = design2)
    return(bcr)
}
#' Plot Bootstrap result
#'
#' Shows the result of a Boostrap with either edge frequencies
#' or confidence intervals
#' @param x bnembs object
#' @param scale numeric value for scaling the edgewidth
#' @param shift numeric value for shifting the edgewidth
#' @param cut shows only edges with a fraction larger than cut
#' @param dec integer for function round
#' @param ci if TRUE shows confidence intervals
#' @param cip range for the confidence interval, e.g. 0.95
#' @param method method to use for conidence itneral
#' computation (see function binom.confint from package binom)
#' @param ... additional parameters for the function mnem::plotDnf
#' from package mnem
#' @author Martin Pirkl
#' @return plot of the network from the bootstrap
#' @method plot bnembs
#' @export
#' @importFrom mnem plotDnf
#' @importFrom binom binom.confint
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t",
#' row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1,
#' maxInhibit = 2, signals = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnemBs(search = "greedy", model = model, CNOlist = CNOlist,
#' fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
#' parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
#' maxSteps = Inf)
plot.bnembs <- function(x, scale = 3, shift = 0.1, cut = 0.5, dec = 2,
                        ci = 0, cip = 0.95, method = "exact", ...) {
    y <- x$x
    n <- x$n
    x <- list()
    x$graph <- sort(unique(y))
    x$freq <- table(y)
    x$freq <- x$freq[order(names(x$freq))]/n
    graph <- x$graph
    freq <- x$freq
    graph <- graph[which(freq >= cut)]
    freq <- freq[which(freq >= cut)]
    freq2 <- NULL
    for (i in seq_len(length(graph))) {
        tmp <- rep(freq[i], length(unlist(strsplit(graph[i], "\\+"))))
        if (length(tmp) > 1) {
            freq2 <- c(freq2, tmp[1])
        }
        freq2 <- c(freq2, tmp)
    }
    freq2 <- as.vector(freq2)
    if (ci) {
        freqn <- freq2*n
        cis <- round(binom::binom.confint(freqn, n, cip,
                                          method = method)[, 5:6], dec)
        cis <-  paste0("[", apply(cis, 1, paste, collapse = " - "), "]")
        plotDnf(graph, edgewidth = freq*scale+shift, edgelabel = cis, ...)
    } else {
        freq2 <- round(freq2, dec)
        plotDnf(graph, edgewidth = freq*scale+shift, edgelabel = freq2, ...)
    }
}
#' Bootstraped Network
#'
#' Runs Bootstraps on the data
#' @param fc matrix with foldchanges or signed effects
#' @param x number of bootstraps
#' @param f percentage to sample, e.g. f = 0.5 samples only 50%
#' the amount of E-genes as the original data
#' @param replace if TRUE classical bootstrap, if FALSE sub-sampling without
#' replacement
#' @param startString matrix with each row being a string denoting a
#' network to start inference several times with a specific network
#' @param ... additional parameters for the bnem function
#' @author Martin Pirkl
#' @return list with the accumulation of edges in x and the number of
#' bootstraps in n
#' @export
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t",
#' row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1,
#' maxInhibit = 2, signals = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnemBs(search = "greedy", model = model, CNOlist = CNOlist,
#' fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
#' parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
#' maxSteps = Inf)
bnemBs <- function(fc, x = 10, f = 0.5, replace = TRUE, startString = NULL,
                   ...) {
    accum <- NULL
    for (i in seq_len(x)) {
        cat(i)
        fcsub <- fc[sample(seq_len(nrow(fc)), ceiling(nrow(fc)*f),
                           replace = replace), ]
        if (is.null(startString)) {
            tmp <- bnem(fc = fcsub, ...)
        } else {
            if (!is(startString, "matrix")) {
                startString <- t(as.matrix(startString))
            }
            score <- 1
            for (j in seq_len(nrow(startString))) {
                tmp0 <- bnem(initBstring = startString[j, ], fc = fcsub, ...)
                if (min(unlist(tmp0$scores)) < score) {
                    tmp <- tmp0
                }
            }
        }

        accum <- c(accum, tmp$graph)
    }

    accum <- list(x = accum, n = x)

    class(accum) <- "bnembs"

    return(accum)
}
#' B-Cell receptor signalling perturbations
#'
#' Processed data from experiments with a stimulated B-Cell receptor (bcr)
#' and perturbed signalling genes. The raw data is available at
#' https://www.ncbi.nlm.nih.gov/geo/ with accession id GSE68761. For the
#' process steps we refer to the publication
#' Martin Pirkl, Elisabeth Hand, Dieter Kube, Rainer Spang, Analyzing
#' synergistic and non-synergistic interactions in signalling pathways
#' using Boolean Nested Effect Models, Bioinformatics, Volume 32, Issue 6,
#' 15 March 2016, Pages 893–900, https://doi.org/10.1093/bioinformatics/btv680.
#' Alternatively see also the function processDataBCR for details.
#' @name bcr
#' @docType data
#' @usage bcr
#' @references Martin Pirkl, Elisabeth Hand, Dieter Kube, Rainer Spang,
#' Analyzing synergistic and non-synergistic interactions in signalling
#' pathways using Boolean Nested Effect Models, Bioinformatics, Volume 32,
#' Issue 6, 15 March 2016, Pages 893–900,
#' https://doi.org/10.1093/bioinformatics/btv680
#' @examples
#' data(bcr)
NA
#' Sample random network and simulate data
#'
#' Draws a random prior network, samples a ground truth from the full boolean
#' extension and generates data
#' @param Sgenes number of S-genes
#' @param maxEdges number of maximum edges in the DAG
#' @param stimGenes number of stimulated S-genes
#' @param layer scaling factor for the sampling of next Sgene layer; high
#' means a higher probability for less Sgenes and low more Sgenes
#' @param frac fraction of hyper-edges in the ground truth
#' @param maxInDeg maximum number of incoming hyper-edges
#' @param dag if TRUE graph will be acyclic
#' @param maxSize maximum number of S-genes in a disjunction or clause
#' @param maxStim maximum of stimulated S-genes in the data samples
#' @param maxInhibit maximum of inhibited S-genes in the data samples
#' @param Egenes number of E-genes per S-gene, e.g. 10 S-genes and 10
#' E-genes will return 100 E-genes overall
#' @param flip number of inhibited E-genes
#' @param reps number of replicates
#' @param keepsif if TRUE does not delete sif file, which encodes
#' the prior network
#' @param maxcount while loopes ensure a reasonable network, maxcount makes sure
#' while loops do not run to infinity
#' @param negation has to be greater or equal to 0 and less than 1; if
#' greater than 0, negation is the sample probability for negative edges
#' @param allstim full network in which all S-genes are possibly stimulated
#' @param and probability for and gates in the GTN
#' @param verbose TRUE for verbose output
#' @author Martin Pirkl
#' @return list with the corresponding prior graph, ground truth network and
#' data
#' @export
#' @importFrom mnem plotDnf
#' @examples
#' sim <- simBoolGtn()
#' plot(sim)
simBoolGtn <-
    function(Sgenes = 10, maxEdges = 25, stimGenes = 2, layer = 1,
             frac = 0.1, maxInDeg = 2,
             dag = TRUE, maxSize = 2, maxStim = 2, maxInhibit = 1,
             Egenes = 10, flip = 0.33, reps = 3, keepsif = FALSE,
             maxcount = 10, negation = 0.25, allstim = FALSE, and = 0.25,
             verbose = FALSE) {
        n <- Sgenes
        m <- Egenes
        p <- layer
        s <- stimGenes
        e <- maxEdges
        r <- reps
        mflip <- flip
        if (allstim) {
            dnf <- randomDnf(n, max.edges = e, max.edge.size = 1, dag = dag,
                             negation = negation)
            cues <- sort(unique(gsub("!", "",
                                     unlist(strsplit(unlist(strsplit(dnf, "=")),
                                                     "\\+")))))
            inputs <-
                unique(unlist(strsplit(gsub("!", "", gsub("=.*", "", dnf)),
                                             "=")))
            outputs <- unique(gsub(".*=", "", dnf))
            stimuli <- cues
            inhibitors <- cues
            both <- stimuli[which(stimuli %in% inhibitors)]
            for (i in both) {
                dnf <- gsub(i, paste(i, "inhibit", sep = ""), dnf)
                dnf <- c(dnf, paste(i, "stim=", i, "inhibit", sep = ""))
                stimuli <- gsub(i, paste(i, "stim", sep = ""), stimuli)
                inhibitors  <- gsub(i, paste(i, "inhibit", sep = ""),
                                    inhibitors)
            }
            pkn <- dnf
        } else {
            Sgenes <- paste0("S", seq_len(n), "g")
            layers <- list()
            layers[[1]] <- stimuli <- prev <- Sgenes[seq_len(s)] # sample(Sgenes, s)
            Sgenes <- inhibitors <- Sgenes[which(!(Sgenes %in% stimuli))]
            pkn <- NULL
            enew <- 0
            count <- 1
            while(length(Sgenes) > 0) {
                count <- count + 1
                pp <- rbeta(1, 1, (length(Sgenes)*p)/10)
                layers[[count]] <- layer <- sample(Sgenes,
                                                   ceiling(length(Sgenes)*pp))
                enew <- enew + length(prev)*length(layer)
                Sgenes <- Sgenes[which(!(Sgenes %in% layer))]
                for (i in seq_len(length(prev))) {
                    if (negation > 0) {
                        die <- sample(c(0,1), 1, prob = c(negation, 1-negation))
                    } else {
                        die <- 1
                    }
                    if (die) {
                        pkn <- c(pkn, paste0(prev[i], "=", layer))
                    } else {
                        pkn <- c(pkn, paste0("!", prev[i], "=", layer))
                    }
                }
                prev <- layer
            }
            dnf <- pkn
            if (enew > e) { e <- enew }
            if (length(dnf) > e) {
                pkn2 <- gsub("!", "", dnf)
                getrid <- sample(which(duplicated(pkn2) == TRUE), length(dnf)-e)
                pkn2 <- dnf[-getrid]
                pkn <- pkn2
                dnf <- pkn
            }
        }
        sifMatrix <- NULL
        for (i in dnf) {
            inputs2 <- unique(unlist(strsplit(gsub("=.*", "", i), "=")))
            output <- unique(gsub(".*=", "", i))
            for (j in inputs2) {
                j2 <- gsub("!", "", j)
                if (j %in% j2) {
                    sifMatrix <- rbind(sifMatrix, c(j, 1, output))
                } else {
                    sifMatrix <- rbind(sifMatrix, c(j2, -1, output))
                }
            }
        }
        mark <- as.numeric(Sys.time())+rnorm(1)
        write.table(sifMatrix, file = paste0("bnemsim", mark, ".sif"),
                    sep = "\t",
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        PKN <- CellNOptR::readSIF(paste0("bnemsim", mark, ".sif"))
        if (!keepsif) {
            unlink(paste0("bnemsim", mark, ".sif"))
        }
        CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors,
                                maxStim = maxStim, maxInhibit = maxInhibit,
                                signals = NULL)
        model <- CellNOptR::preprocessing(CNOlist, PKN,
                                          maxInputsPerGate=maxSize,
                                          verbose = verbose)
        bString <- numeric(length(model$reacID))
        if (allstim) {
            bString[sample(seq_len(length(model$reacID)),
                           ceiling(length(model$reacID)*frac))] <- 1
        } else {
            for (i in seq_len(length(layers))) {
                for (j in seq_len(length(layers[[i]]))) {
                    if (i == 1) {
                        gates <- grep(paste0(layers[[i]][j]),
                                      model$reacID)
                    } else {
                        gates <- grep(paste0("=", layers[[i]][j], "$"),
                                      model$reacID)
                    }
                    if (length(gates) == 1) {
                        bString[gates] <- 1
                    } else {
                        prob <- rep(1-and, length(gates))
                        prob[grep("\\+", model$reacID[gates])] <- and
                        bString[sample(gates,
                                       max(c(floor(length(gates)*frac), 1)),
                                       prob = prob)] <- 1
                    }
                }
            }
            toomany <- table(gsub(".*=", "", model$reacID[as.logical(bString)]))
            toomany <- toomany[which(toomany > maxInDeg)]
            for (i in seq_len(length(toomany))) {
                manyIn <- intersect(which(bString == 1), grep(paste0("=", names(toomany)[i]), model$reacID))
                bString[grep(paste0("=", names(toomany)[i]), model$reacID)] <- 0
                bString[sample(manyIn, maxInDeg)] <- 1
            }
        }
        bString <- reduceGraph(bString, model, CNOlist)
        steadyState <- steadyState2 <- simulateStatesRecursive(CNOlist, model,
                                                               bString)
        ind <- grep(paste(inhibitors, collapse = "|"), colnames(steadyState2))
        steadyState2[, ind] <- steadyState2[, ind] + CNOlist@inhibitors
        exprs <- t(steadyState)[rep(seq_len(ncol(steadyState)), m),
                                rep(seq_len(nrow(steadyState)), r)]
        ERS <- computeFc(CNOlist, t(steadyState))
        stimcomb <- apply(expand.grid(stimuli, stimuli), c(1,2), as.character)
        stimuli.pairs <- apply(stimcomb, 1, paste, collapse = "_")
        ind <- grep(paste(c(paste("Ctrl_vs_", c(stimuli, inhibitors), sep = ""),
                            paste(stimuli, "_vs_", stimuli, "_",
                                  rep(inhibitors, each = length(stimuli)),
                                  sep = ""),
                            paste(stimuli.pairs, "_vs_", stimuli.pairs, "_",
                                  rep(inhibitors, each = length(stimuli.pairs)),
                                  sep = "")),
                          collapse = "|"), colnames(ERS))
        ERS <- ERS[, ind]
        fc <- ERS[rep(seq_len(nrow(ERS)), m), rep(seq_len(ncol(ERS)), r)]
        flip <- sample(seq_len(nrow(fc)), floor(mflip*row(fc)))
        fc[flip, ] <- fc[flip, ]*(-1)
        rownames(fc) <- paste(rownames(fc), seq_len(nrow(fc)), sep = "_")
        sim <- list(PKN = PKN, CNOlist = CNOlist, model = model,
                    bString = bString, fc = fc, exprs = exprs, ERS = ERS)
        class(sim) <- "bnemsim"
        return(sim)
    }
#' Inverse absorption
#'
#' applies "inverse" absorption law to a disjuncitve normal form
#' @param bString a disjunctive normal form or binary vector according to model
#' @param model model for respective binary vector
#' @author Martin Pirkl
#' @return bString after "inverse" absorption law
#' @export
#' @import mnem
#' @examples
#' graph <- c("A+B=C", "A=C")
#' absorptionII(graph)
absorptionII <-
    function(bString, model=NULL) {
        if (is.null(model)) {
            graph <- bString
            nodes <-
                unique(gsub("!", "",
                            unlist(strsplit(unlist(strsplit(graph, "=")),
                                            "\\+"))))
        } else {
            graph <- model$reacID[which(bString == 1)]
            nodes <- model$namesSpecies
        }
        for (i in graph) {
            players <- unlist(strsplit(gsub("=.*", "", i), "\\+"))
            target <- gsub(".*=", "", i)
            others <- nodes[-which(nodes %in% c(players, target))]
            players2 <- gsub("!", "", players)
            change1 <- which(players == players2)
            change2 <- which(!(players == players2))
            if (length(change1) > 0) {
                others <- c(others, paste("!", players2[change1],
                                          sep = ""))
            }
            if (length(change2) > 0) {
                others <- c(others[-which(others %in% players2[change2])],
                            paste("\\+", players2[change2], sep = ""),
                            paste("^", players2[change2], sep = ""))
            }
            targets <-
                intersect(grep(paste(paste(
                    "^", paste(players, collapse = "|^"), sep = ""), "|",
                    paste("+", paste(players, collapse = "|+"), sep = ""),
                    sep = ""), graph), grep(paste("=", target, sep = ""),
                                            graph))
            toomuch <- which(targets %in% grep(paste(others, collapse = "|"),
                                               graph))
            if (length(toomuch) > 0) {
                targets <- targets[-toomuch]
            }
            if (length(targets) > 1) {
                targets <- targets[-which(targets %in% which(graph %in% i))]
                if (is.null(model)) {
                    if (sum(bString %in% graph[targets]) > 0) {
                        bString <- bString[-which(bString %in% graph[targets])]
                    }
                } else {
                    bString[which(model$reacID %in% graph[targets])] <- 0
                }
            }
        }
        return(bString)
    }
#' Absorption
#'
#' applies absorption law to a disjuncitve normal form
#' @param bString a disjunctive normal form or binary vector according to model
#' @param model model for respective binary vector
#' @author Martin Pirkl
#' @return bString after absorption law
#' @export
#' @examples
#' graph <- c("A+B=C", "A=C")
#' absorption(graph)
absorption <-
    function(bString, model=NULL) {
        if (is.null(model)) {
            graph <- bString
        } else {
            graph <- model$reacID[which(bString == 1)]
        }
        for (i in graph) {
            targets <-
                grep(paste("(?=.*", gsub("\\+", ")(?=.*",
                                         gsub("=", ")(?=.*=", i)), ")",
                           sep = ""), graph, perl = TRUE)
            toomuch <-
                grep(paste("!", gsub("\\+", "|!",
                                     gsub("=.*", "", i)), "",
                           sep = ""), graph[targets])
            if (length(toomuch) > 0) {
                targets <-
                    targets[-grep(paste("!", gsub("\\+", "|!",
                                                  gsub("=.*", "", i)), "",
                                        sep = ""), graph[targets])]
            }
            if (length(targets) > 1) {
                targets <- targets[-which(targets == which(graph %in% i))]
                if (is.null(model)) {
                    if (sum(bString %in% graph[targets]) > 0) {
                        bString <- bString[-which(bString %in% graph[targets])]
                    }
                } else {
                    bString[which(model$reacID %in% graph[targets])] <- 0
                }
            }
        }
        return(bString)
    }
#' Boolean Nested Effects Model main function
#'
#' This function takes a prior network and normalized perturbation data
#' as input and trains logical functions on that prior network
#' @param search Type of search heuristic. Either "greedy", "genetic" or
#' "exhaustive". "greedy" uses a greedy algorithm to move through the local
#' neighbourhood of a initial hyper-graph. "genetic" uses a genetic algorithm.
#' "exhaustive" searches through the complete search space and is not
#' recommended.
#' @param fc Foldchanges of gene expression values or equivalent input
#' (normalized pvalues, logodds, ...). If left NULL, the gene expression
#' data is used to calculate naive foldchanges.
#' @param exprs Optional normalized gene expression data.
#' @param egenes list object; each list entry is named after an S-gene and
#' contains the egenes which are potential children
#' @param pkn Prior knowledge network.
#' @param design Optional design matrix for the gene expression values, if
#' available. If kept NULL, bnem needs either stimuli, inhibitors or a CNOlist
#' object.
#' @param stimuli Character vector of stimuli names.
#' @param inhibitors Character vector of inhibitors.
#' @param signals Optional character vector of signals. Signals are S-genes,
#' which can directly regulate E-genes. If left NULL, alls stimuli and
#' inhibitors are defined as signals.
#' @param CNOlist CNOlist object if available.
#' @param model Model object including the search space, if available.
#' @param sizeFac Size factor penelizing the hyper-graph size.
#' @param NAFac factor penelizing NAs in the data.
#' @param parameters parameters for discrete case (not recommended);
#' has to ba list with entries cutOffs and scoring:
#' cutOffs = c(a,b,c) with a (cutoff for real zeros),
#' b (cutoff for real effects),
#' c = -1 for normal scoring, c between 0 and
#' 1 for keeping only relevant % of E-genes,
#' between -1 and 0 for keeping only a specific quantile of E-genes,
#' and c > 1 for keeping the top c E-genes;
#' scoring = c(a,b,c) with a (weight for real effects),
#' c (weight for real zeros),
#' b (multiplicator for effects/zeros between a and c);
#' @param parallel Parallelize the search. An integer value specifies the
#' number of threads on the local machine. A list object as in list(c(1,2,3),
#' c("machine1", "machine2", "machine3")) specifies the threads distributed
#' on different machines (local or others).
#' @param method Scoring method can be "cosine" a correlation,
#' distance measure or a
#' probability based score "llr". See ?cor and ?dist for details.
#' @param relFit if TRUE a relative fit for each
#' E-gene is computed (not recommended)
#' @param verbose TRUE gives additional information during the search.
#' @param reduce if TRUE reduces the search space for exhaustive search
#' @param initBstring Binary string of the initial hyper-graph.
#' @param popSize Population size (only "genetic").
#' @param pMutation Probability for mutation (only "genetic").
#' @param maxTime Define a maximal time (seconds) for the search.
#' @param maxGens Maximal number of generations (only "genetic").
#' @param stallGenMax Maximum number of stall generations (only "genetic").
#' @param relTol Score tolerance for networks defined as optimal but with a
#' lower score as the real optimum (only "genetic").
#' @param priorBitString Binary string defining hyper-edges which are added
#' to every hyper-graph. E.g. if you know hyper-edge 55 is definitly there and
#' to fix that, set priorBitString[55] <- 1 (only "genetic").
#' @param selPress Selection pressure for the stochastic universal sampling
#' (only "genetic").
#' @param approach default "fc" for foldchanges or signed effect probabilities,
#' "abs" for absolute effects or probabilities
#' @param fit "linear" or "nonlinear fit for stochastic universal sampling
#' @param targetBstring define a binary string representing a network;
#' if this network is found, the computation stops
#' @param elitism Number of best hyper-graphs transferred to the next generation
#' (only "genetic").
#' @param inversion Number of worst hyper-graphs for which their binary strings
#' are inversed  (only "genetic").
#' @param parallel2 if TRUE parallelises the starts and not the search itself
#' @param selection "t" for tournament selection and "s" for stochastic
#' universal sampling (only "genetic").
#' @param type type of the paralellisation on multpile machines (default: SOCK)
#' @param exhaustive If TRUE an exhaustive search is conducted if the genetic
#' algorithm would take longer (only "genetic").
#' @param delcyc If TRUE deletes cycles in all hyper-graphs (not recommended).
#' @param seeds how many starts for the greedy search? (default: 1)
#' @param maxSteps Maximal number of steps (only "greedy").
#' @param node vector of S-gene names, which are used in the greedy
#' search; if node = NULL all nodes are considered
#' @param absorpII Use inverse absorption.
#' @param draw If TRUE draws the network evolution.
#' @param prior Binary vector. A 1 specifies hyper-edges which should not be
#' optimized (only "greedy").
#' @param maxInputsPerGate If no model is supplied, one is created with
#' maxInputsPerGate as maximum number of parents for each hyper-edge.
#' @author Martin Pirkl
#' @seealso nem
#' @export
#' @import
#' CellNOptR
#' snowfall
#' mnem
#' methods
#' @return List object including the optimized hyper-graph and its
#' corresponding binary string.
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t",
#' row.names = FALSE, col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1,
#' maxInhibit = 2, signals = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnem(search = "greedy", model = model, CNOlist = CNOlist,
#' fc = fc, pkn = PKN, stimuli = "A", inhibitors = c("B","C","D"),
#' parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
#' maxSteps = Inf)
bnem <-
    function(search = "greedy",

             fc=NULL,
             exprs=NULL,
             egenes=NULL,
             pkn=NULL,
             design=NULL,
             stimuli=NULL,
             inhibitors=NULL,
             signals=NULL,
             CNOlist=NULL,
             model=NULL,
             sizeFac=10^-10,
             NAFac=1,
             parameters=list(cutOffs=c(0,1,0), scoring=c(0.1,0.2,0.9)),
             parallel = NULL,
             method = "s",
             approach = "fc",
             relFit = FALSE,
             verbose = TRUE,
             reduce = TRUE,
             parallel2 = 1,

             initBstring = NULL,
             popSize = 100,
             pMutation = 0.5,
             maxTime = Inf,
             maxGens = Inf,
             stallGenMax = 10,
             relTol = 0.01,
             priorBitString = NULL,
             selPress = c(1.2,0.0001),
             fit = "linear",
             targetBstring = "none",
             elitism = NULL,
             inversion = NULL,
             selection = c("t"),
             type = "SOCK",
             exhaustive = FALSE,
             delcyc = FALSE,

             seeds = 1,
             maxSteps = Inf,
             node = NULL,
             absorpII = TRUE,
             draw = TRUE,
             prior = NULL,
             maxInputsPerGate = 2
             ) {
        if (is.null(model) | is.null(CNOlist)) {
            tmp <- preprocessInput(stimuli=stimuli,inhibitors=inhibitors,
                                   signals=signals,design=design,exprs=exprs,
                                   fc=fc,pkn=pkn,
                                   maxInputsPerGate=maxInputsPerGate)

            CNOlist <- tmp$CNOlist
            NEMlist <- tmp$NEMlist
            model <- tmp$model
        } else {
            NEMlist <- list()
            NEMlist$fc <- fc
            NEMlist$exprs <- exprs
            NEMlist$egenes <- egenes
        }
        CNOlist <- checkCNOlist(CNOlist)

        if (search %in% c("greedy", "genetic", "exhaustive")) {

            if (search %in% "greedy") {
                res <- localSearch(CNOlist=CNOlist, NEMlist=NEMlist,
                                   model=model,
                                   approach = approach, initSeed = initBstring,
                                   seeds = seeds, parameters = parameters,
                                   sizeFac = sizeFac, NAFac = NAFac,
                                   relTol = relTol, verbose = verbose,
                                   parallel=parallel, parallel2 = parallel2,
                                   relFit = relFit, method = method,
                                   max.steps = maxSteps, max.time = maxTime,
                                   node = node, absorpII = absorpII,
                                   draw = draw,
                                   prior = prior)
                minSeed <- 1
                if (length(res$scores) > 1) {
                    for (i in seq_len(length(res$scores))) {
                        if (min(res$scores[[i]]) < min(res$scores[[minSeed]])) {
                            minSeed <- i
                        }
                    }
                }
                bString <- res$bStrings[minSeed, ]
                result <- list(graph = model$reacID[as.logical(bString)],
                               bString = bString, bStrings = res$bStrings,
                               scores = res$scores)
            }
            if (search %in% "genetic") {
                res <- gaBinaryNemT1(CNOlist=CNOlist, model=model,
                                     initBstring = initBstring,
                                     sizeFac = sizeFac,
                                     NAFac = NAFac,popSize = popSize,
                                     pMutation = pMutation,maxTime = maxTime,
                                     maxGens = maxGens,
                                     stallGenMax = stallGenMax,relTol = relTol,
                                     verbose = verbose,
                                     priorBitString = priorBitString,
                                     selPress = selPress,approach = approach,
                                     NEMlist=NEMlist,fit = fit,
                                     targetBstring = targetBstring,
                                     elitism = elitism,inversion = inversion,
                                     graph = draw,parameters = parameters,
                                     parallel = parallel,parallel2 = parallel2,
                                     selection = selection,relFit = relFit,
                                     method = method,type = type,
                                     exhaustive = exhaustive,delcyc = delcyc
                                     )
                result <- list(graph = model$reacID[as.logical(res$bString)],
                               bString = res$bString, bStrings = res$stringsTol,
                               scores = res$stringsTolScores)
            }
            if (search %in% "exhaustive") {
                res <- exSearch(CNOlist=CNOlist,model=model,sizeFac=sizeFac,
                                NAFac=NAFac,NEMlist=NEMlist,
                                parameters=parameters, parallel = parallel,
                                method = method, relFit = relFit,
                                verbose = verbose, reduce = reduce,
                                approach = approach)
                result <- list(graph = model$reacID[as.logical(res$bString)],
                               bString = res$bString, bStrings = res$bStrings,
                               scores = res$scores)
            }
            return(result)
        } else {
            return("search must be be one of greedy, genetic or exhaustive")
        }
    }
#' Compute repsonse scheme
#'
#' computes response scheme given an activation pattern
#' (absolute gene expression, truth table)
#' @param CNOlist a CNOlist object with correct annotation
#' @param y activation pattern according to the annotation in CNOlist
#' @author Martin Pirkl
#' @return numeric matrix with annotated response scheme
#' @export
#' @import CellNOptR
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signals = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
computeFc <-
    function (CNOlist, y) {
        test <- 1 # for debugging
        CompMat <- numeric()
        CompMatNames <- character()
        cnolistStimuli <- apply(CNOlist@stimuli, 1, sum)
        cnolistInhibitors <- apply(CNOlist@inhibitors, 1, sum)
        cnolistCues <- apply(CNOlist@cues, 1, sum)
        maxStim <- max(cnolistStimuli)
        maxKd <- max(cnolistInhibitors)
        grepCtrl <- which(cnolistCues == 0)[1]
        grepStims <- intersect(which(cnolistStimuli != 0),
                               which(cnolistInhibitors == 0))
        grepKds <- intersect(which(cnolistStimuli == 0),
                             which(cnolistInhibitors != 0))
        grepStimsKds <- intersect(which(cnolistStimuli != 0),
                                  which(cnolistInhibitors != 0))
        stimsKdsCbind <- cbind(CNOlist@stimuli, CNOlist@inhibitors)
        ## get ctrl_vs_kd:
        inhibitorsNames <- NULL
        for (i in grepKds) {
            inhibitorsNames <-
                c(inhibitorsNames,
                  paste(colnames(CNOlist@inhibitors)[which(
                                    CNOlist@inhibitors[i, ] >= 1)],
                        collapse = "_"))
        }
        if (length(grepKds) > 0) {
            CompMat <- cbind(CompMat, y[, grepKds] - y[, grepCtrl])
            CompMatNames <-
                c(CompMatNames, paste("Ctrl_vs_", inhibitorsNames,
                                      sep = ""))
        }
        ## get ctrl_vs_stim:
        stimuliNames <- NULL
        for (i in grepStims) {
            if (sum(CNOlist@stimuli[i, ] >= 1) > 1) {
                stimuliNames <-
                    c(stimuliNames,
                      paste(names(which(CNOlist@stimuli[i, ] >= 1)),
                            collapse = "_"))
            } else {
                stimuliNames <-
                    c(stimuliNames,
                      colnames(CNOlist@stimuli)[
                          which(CNOlist@stimuli[i, ] >= 1)])
            }
        }
        if (length(grepStims) > 0) {
            CompMat <- cbind(CompMat, y[, grepStims] - y[, grepCtrl])
            CompMatNames <- c(CompMatNames, paste("Ctrl_vs_", stimuliNames,
                                                  sep = ""))
        }
        ## get stim_vs_stim_kd:
        ## combiNames <- NULL
        ## for (i in grepStimsKds) {
        ##     combiNames <-
        ##         c(combiNames,
        ##           paste(names(which(stimsKdsCbind[i, ] >= 1)),
        ##                 collapse = "_"))
        ## }
        combiNames <- rownames(CNOlist@cues)[grepStimsKds]
        if (length(grepStimsKds) > 0 & length(grepStims) > 0) {
            CompMat <-
                cbind(CompMat,
                      y[, rep(grepStimsKds, length(grepStims))] -
                      y[, sort(rep(grepStims, length(grepStimsKds)))])
            orderStims <- order(rep(grepStims, length(grepStimsKds)))
            CompMatNames <-
                c(CompMatNames,
                  paste(rep(stimuliNames, length(combiNames))[orderStims],
                        "_vs_", rep(combiNames, length(stimuliNames)),
                        sep = ""))
        }
        ## combine:
        colnames(CompMat) <- CompMatNames
        if (sum(duplicated(colnames(CompMat)) == TRUE)) {
            CompMat <-
                CompMat[, -which(duplicated(colnames(CompMat)) == TRUE)]
        }
        if (test == 1) {
            ## get stim_vs_stim:
            combiNames2 <- NULL
            for (i in grepStims) {
                combiNames2 <-
                    c(combiNames2,
                      paste(names(which(CNOlist@stimuli[i, ] >= 1)),
                            collapse = "_"))
            }
            if (length(grepStims) > 0) {
                CompMat <-
                    cbind(CompMat, y[, rep(grepStims, length(grepStims))] -
                                   y[, sort(rep(grepStims, length(grepStims)))])
                orderStims2 <- order(rep(grepStims, length(grepStims)))
                CompMatNames <-
                    c(CompMatNames,
                      paste(rep(stimuliNames,
                                length(combiNames2))[orderStims2], "_vs_",
                            rep(combiNames2, length(stimuliNames)), sep = ""))
            }
            ## get kd_vs_stim_kd:
            combiNames <- NULL
            for (i in grepStimsKds) {
                combiNames <-
                    c(combiNames,
                      paste(names(which(stimsKdsCbind[i, ] >= 1)),
                            collapse = "_"))
            }
            if (length(grepStimsKds) > 0 & length(grepKds) > 0) {
                CompMat <-
                    cbind(CompMat, y[, rep(grepStimsKds, length(grepKds))] -
                                   y[, sort(rep(grepKds,
                                                length(grepStimsKds)))])
                orderKds <- order(rep(grepKds, length(grepStimsKds)))
                CompMatNames <-
                    c(CompMatNames,
                      paste(rep(inhibitorsNames,
                                length(combiNames))[orderKds], "_vs_",
                            rep(combiNames, length(inhibitorsNames)), sep = ""))
            }
            ## get ctrl_vs_stim_kd:
            combiNames <- NULL
            for (i in grepStimsKds) {
                combiNames <-
                    c(combiNames,
                      paste(names(which(stimsKdsCbind[i, ] >= 1)),
                            collapse = "_"))
            }
            if (length(grepStimsKds) > 0) {
                CompMat <- cbind(CompMat, y[, grepStimsKds] - y[, grepCtrl])
                CompMatNames <- c(CompMatNames, paste("Ctrl_vs_", combiNames,
                                                      sep = ""))
            }
        }
        ## combine:
        colnames(CompMat) <- CompMatNames
        if (sum(duplicated(colnames(CompMat)) == TRUE)) {
            CompMat <-
                CompMat[, -which(duplicated(colnames(CompMat)) == TRUE)]
        }
        return(CompMat)
    }
#' Convert normal form
#'
#' converts a disjunctive normal form into a conjunctive normal form and
#' vice versa;
#' input graph as disjunctive normal form like that:
#' c("A+B=D", "C=D", "G+F=U", ...); output is the dual element
#' also in disjunctive normal form;
#' @param g graph in normal form
#' @author Martin Pirkl
#' @return converted graph normal form
#' @export
#' @examples
#' g <- "A+B=C"
#' g2 <- convertGraph(g)
convertGraph <-
    function(g) {
        g <- sort(g)
        targets <- gsub(".*=", "", g)
        g.new <- NULL
        for (i in unique(targets)) {
            dnf <- list()
            count <- 1
            for (j in g[grep(paste("=", i, sep = ""), g)]) {
                dnf[[count]] <-
                    sort(unique(unlist(strsplit(gsub("=.*", "", j), "\\+"))))
                count <- count + 1
            }
            cnf <- expand.grid(dnf)
            dnf <- NULL
            for (j in seq_len(dim(cnf)[1])) {
                dnf <- c(dnf, paste(sort(unique(unlist(cnf[j, ]))),
                                    collapse = "+"))
            }
            dnf <- paste(sort(dnf), "=", i, sep = "")
            g.new <- c(g.new, dnf)
        }
        vertices <- sort(unique(unlist(strsplit(unlist(strsplit(g.new, "=")),
                                                "\\+"))))
        for (i in vertices) {
            if (length(grep(paste(i, ".*", i, ".*=", sep = ""), g.new)) > 0) {
                g.new <- g.new[-grep(paste(i, ".*", i, ".*=", sep = ""), g.new)]
            }
        }
        return(g.new)
    }
#' Create dummy CNOlist
#'
#' creates a general CNOlist object from meta information
#' @param stimuli character vector of stimulated genes
#' @param inhibitors character vector of inhibited genes
#' @param maxStim maximal number of stimulated genes for a single experiment
#' @param maxInhibit maximal number of inhibited genes for a single experiment
#' @param signals character vector of genes which can directly regulate effect
#' reporters
#' @author Martin Pirkl
#' @return CNOlist object
#' @export
#' @import CellNOptR
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signals = c("A", "B","C","D"))
dummyCNOlist <-
    function(stimuli = NULL, inhibitors = NULL, maxStim = 0, maxInhibit = 0,
             signals = NULL) {
        if (is.null(signals)) {
            signals <- c(stimuli, inhibitors)
        }
        stimn <- length(stimuli)
        inhibn <- length(inhibitors)
        ## do the stim design:
        if (maxStim > stimn) {
            maxStim <- stimn
        }
        if (stimn > 0 & maxStim > 0) {
            mat.size <- 0
            for (i in seq_len(maxStim)) {
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
                    combs <- t(combn(seq_len(stimn), i))
                    for (j in seq_len(nrow(combs))) {
                        count <- count + 1
                        stimDesign[count, combs[j, ]] <- 1
                    }
                }
            }
            colnames(stimDesign) <- stimuli
        }
        ## do the inhib design:
        inhibn <- length(inhibitors)
        if (maxInhibit > inhibn) {
            maxInhibit <- inhibn
        }
        if (inhibn > 0 & maxInhibit > 0) {
            mat.size <- 0
            for (i in seq_len(maxInhibit)) {
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
                    combs <- t(combn(seq_len(inhibn), i))
                    for (j in seq_len(nrow(combs))) {
                        count <- count + 1
                        inhibDesign[count, combs[j, ]] <- 1
                    }
                }
            }
            colnames(inhibDesign) <- inhibitors
        }
        ## put them together in combinations:
        if (stimn > 0 & inhibn > 0) {
            if (maxStim > 0 & maxInhibit > 0) {
                design <- matrix(0, nrow(stimDesign)*nrow(inhibDesign),
                                 stimn+inhibn)
                for (i in seq_len(nrow(stimDesign))) {
                    design[((i-1)*nrow(inhibDesign) + 1):
                           (i*nrow(inhibDesign)), ] <-
                        cbind(stimDesign[rep(i, nrow(inhibDesign)), , drop = FALSE],
                              inhibDesign)
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
        (ncol(design) - ncol(stimDesign)))), cbind(matrix(0, nrow(inhibDesign),
        (ncol(design) - ncol(inhibDesign))), inhibDesign), design)
        colnames(design) <- colnamesdesign
        ## make signalmatrix:
        signaln <- length(signals)
        if (signaln > 0) {
            signalData <- matrix(0, nrow(design)+1, signaln)
            colnames(signalData) <- signals
        } else {
            signalData <- matrix(0, nrow(design)+1, 1)
        }
        smult <- nrow(design)/nrow(stimDesign)
        imult <- nrow(design)/nrow(inhibDesign)
        design <- rbind(0, design)
        stimDesign <- as.matrix(design[, seq_len(ncol(stimDesign))])
        inhibDesign <- as.matrix(design[, (ncol(stimDesign)+1):ncol(design)])
        rownames(design) <- rownames(inhibDesign) <- rownames(stimDesign) <-
            rownames(signalData) <- c("Ctrl", 2:nrow(design))
        getRowname <- function(i, M) {
            r <- paste(colnames(M)[which(M[i, ] == 1)], collapse = "_")
            return(r)
        }
        rownames(design)[2:nrow(design)] <-
            rownames(inhibDesign)[2:nrow(design)] <-
            rownames(stimDesign)[2:nrow(design)] <-
            rownames(signalData)[2:nrow(design)] <-
            unlist(lapply(as.list(2:nrow(design)), getRowname, design))
        if (ncol(stimDesign) == 1) {
            colnames(stimDesign) <- stimuli
        }
        if (ncol(inhibDesign) == 1) {
            colnames(inhibDesign) <- inhibitors
        }
        cnolist <- new("CNOlist",
                       cues = design, inhibitors = inhibDesign,
                       stimuli = stimDesign,
                       signals = list(signalData, signalData),
                       timepoints = as.character(c(0,1)))
        cnolist <- checkCNOlist(cnolist)
        return(cnolist)
    }
#' Switch between epiNEM and B-NEM
#'
#' Convert epiNEM model into general Boolean graph.
#' Only needed for comparing accuracy of inferred network for bnem and epiNEM.
#' @param t full epiNEM model
#' @author Martin Pirkl
#' @seealso CreateTopology
#' @export
#' @import epiNEM
#' @examples
#' topology <- epiNEM::CreateTopology(3, 1, force = TRUE)
#' topology <- unlist(unique(topology), recursive = FALSE)
#' extTopology <- epiNEM::ExtendTopology(topology$model, 100)
#' b <- epiNEM2Bg(extTopology)
#' @return differential effects pattern
epiNEM2Bg <- function(t) {
    if (is.matrix(t)) {
        colnames(t) <- paste("S_vs_S_", gsub("\\.", "_", colnames(t)),
                             sep = "")
        return(t)
    } else {
        tmp <- apply(t$origModel, 2, sum)
        stim <- rownames(t$origModel)[which(tmp == min(tmp))]
        graph <- NULL

        for (i in seq_len(length(t$column))) {
            parents <-
                sort(rownames(t$origModel)[which(t$origModel[
                                                     , t$column[i]] == 1)])
            child <- colnames(t$origModel)[t$column[i]]
            if (length(parents) == 2) {
                if (t$logics[i] %in% "OR") {
                    graph <- unique(c(graph,
                                      convertGraph(adj2dnf(t$origModel))))
                }
                if (t$logics[i] %in% "AND") {
                    graph <-
                        unique(c(graph,
                                 transRed(
                                     convertGraph(adj2dnf(t$origModel)))))
                    graph <-
                        c(graph,
                          convertGraph(graph[
                              grep(paste(paste(parents,
                                               collapse = ".*\\+.*"),
                                         child, sep = "="), graph)]))
                    graph <-
                        graph[-grep(paste(paste(parents,
                                                collapse = ".*\\+.*"),
                                          child, sep = "="), graph)]
                }
                if (t$logics[i] %in% paste(parents[2], " masks the effect of ",
                                           parents[1], sep = "")) {
                    graph <- c(graph,
                               unique(convertGraph(adj2dnf(t$origModel))))
                    graph <-
                        c(graph,
                          gsub(parents[2],
                               paste("!", parents[2], sep = ""),
                               gsub("\\+\\+|^\\+", "",
                                    gsub(parents[1], "",
                                         graph[
                                             grep(paste(paste(
                                                 parents,
                                                 collapse =
                                                     ".*\\+.*"), child,
                                                 sep = "="), graph)]))),
                          gsub("\\+=", "=",
                               gsub("\\+\\+|^\\+", "",
                                    gsub(parents[2], "",
                                         graph[grep(paste(paste(parents,
                                                                collapse =
                                                                    ".*\\+.*"),
                                                          child, sep = "="),
                                                    graph)]))))
                    graph <- graph[-grep(paste(paste(parents,
                                                     collapse = ".*\\+.*"),
                                               child, sep = "="), graph)]
                }
                if (t$logics[i] %in% paste(parents[1],
                                           " masks the effect of ",
                                           parents[2], sep = "")) {
                    graph <- c(graph,
                               unique(convertGraph(adj2dnf(t$origModel))))
                    graph <-
                        c(graph,
                          gsub(parents[1],
                               paste("!", parents[1], sep = ""),
                               gsub("\\+\\+|^\\+", "",
                                    gsub(parents[2], "",
                                         graph[grep(paste(paste(parents,
                                                                collapse =
                                                                    ".*\\+.*"),
                                                          child, sep = "="),
                                                    graph)]))),
                          gsub("\\+=", "=",
                               gsub("\\+\\+|^\\+", "",
                                    gsub(parents[1], "",
                                         graph[grep(paste(paste(parents,
                                                                collapse =
                                                                    ".*\\+.*"),
                                                          child, sep = "="),
                                                    graph)]))))
                    graph <- gsub("\\+=", "=", graph)
                    graph <-
                        graph[-grep(paste(paste(parents,
                                                collapse = ".*\\+.*"),
                                          child, sep = "="), graph)]
                }
                if (t$logics[i] %in% "XOR") {
                    graph <- unique(c(graph,
                                      convertGraph(adj2dnf(t$origModel))))
                    edge <-  graph[grep(paste(paste(parents,
                                                    collapse = ".*\\+.*"),
                                              child, sep = "="), graph)]
                    for (j in parents) {
                        edge <- gsub(j, paste("!", j, sep = ""), edge)
                    }
                    graph <- unique(c(graph, edge))
                }
            }
            if (length(parents) > 2 | length(parents) == 1) {
                graph <- c(graph, paste(sort(parents), child, sep = "="))
            }
        }
        all <- rownames(t$origModel)
        children2 <- unique(gsub(".*=", "", graph))
        if (sum(!(all %in% children2)) > 0) {
            graph <- c(graph, paste("S", all[which(!(all %in% children2))],
                                    sep = "="))
        }
        return(unique(graph))
    }
}
#' compute residuals
#'
#' calculates residuals (data and optimized network do not match) and
#' visualizes them
#' @param bString Binary vector denoting the network given a model
#' @param CNOlist CNOlist object
#' @param model Network model object.
#' @param fc ORS of the data as numeric matrix.
#' @param exprs Optional activation scheme of the data as numeric matrix.
#' @param egenes Atachment of the E-genes (optional) as
#' list named after S-genes.
#' @param NEMlist NEMlist object (optional).
#' @param parameters see ?bnem
#' @param method Scoring method (optional).
#' @param sizeFac Zeta parameter to penelize network size.
#' @param main Main title of the figure.
#' @param sub Subtitle of the figure.
#' @param cut If TRUE does not visualize experiments/S-genes which do
#' not have any residuals.
#' @param approach see ?bnem
#' @param parallel the number of threads used for computation.
#' @param verbose verbose output
#' @param ... additional parameters for ?epiNEM::HeatmapOP
#' @author Martin Pirkl
#' @return numeric matrices indicating experiments and/or genes, where the
#' network and the data disagree
#' @export
#' @import
#' CellNOptR
#' snowfall
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signal = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnem(search = "greedy", CNOlist = CNOlist, fc = fc, model = model,
#' parallel = NULL, initBstring = initBstring, draw = FALSE, verbose = FALSE,
#' maxSteps = Inf)
#' rownames(fc) <- seq_len(nrow(fc))
#' ## val <- validateGraph(CNOlist = CNOlist, fc = fc, model = model,
#' ## bString = res$bString, Egenes = 10, Sgene = 4)
#' residuals <- findResiduals(res$bString, CNOlist, model, fc = fc)
findResiduals <-
    function(bString, CNOlist, model, fc=NULL, exprs=NULL, egenes=NULL,
             NEMlist=NULL, parameters = list(cutOffs = c(0,1,0),
                                             scoring = c(0.1,0.2,0.9)),
             method = "s", sizeFac = 10^-10,
             main = "residuals for decoupled vertices",
             sub = paste0("green residuals are added effects (left positive,",
                          " right negative) and red residuals are deleted ",
                          "effects"),
cut = TRUE, approach = "fc", parallel = NULL, verbose = TRUE, ...) {
        if (is.null(NEMlist)) {
            NEMlist <- list()
            NEMlist$fc <- fc
            NEMlist$egenes <- egenes
            NEMlist$exprs <- exprs
        }
        CNOlist <- checkCNOlist(CNOlist)
        method <- checkMethod(method)[1]
        NEMlist <- checkNEMlist(NEMlist, CNOlist, parameters, approach, method)
        simResults <- simulateStatesRecursive(CNOlist = CNOlist, model = model,
                                              bString = bString)
        SCompMat <- computeFc(CNOlist, t(simResults))
        SCompMat <- SCompMat[, which(colnames(SCompMat) %in%
                                     colnames(NEMlist$fc))]
        NEMlist$fc <- NEMlist$fc[, order(colnames(NEMlist$fc))]
        SCompMat <- SCompMat[, order(colnames(SCompMat))]
        SCompMat <- SCompMat[, colnames(NEMlist$fc)]
        stimuli <- colnames(CNOlist@stimuli)
        inhibitors <- colnames(CNOlist@inhibitors)
        tmp <- computeScoreNemT1(CNOlist, model = model, bString,
                                 NEMlist = NEMlist, tellme = 1,
                                 parameters = parameters, method = method,
                                 verbose = verbose, sizeFac = sizeFac)
        EtoS <- tmp$EtoS

        if (verbose) {
            print(paste("calculating residuals for ",
                        ncol(CNOlist@signals[[1]]),
                        " S-genes based on ", length(unique(EtoS[, 1])),
                        " E-genes.", sep = ""))
        }

        resMat <- matrix(0, nrow = ncol(CNOlist@signals[[1]]),
                         ncol = 2*ncol(NEMlist$fc))
        resVec <- numeric(ncol(CNOlist@signals[[1]]))
        resType <- matrix(0, nrow = ncol(CNOlist@signals[[1]]),
                          ncol = 2*ncol(NEMlist$fc))

        checkSgene <- function(i) {
            resType <- numeric(2*ncol(NEMlist$fc))
            resMat <- numeric(2*ncol(NEMlist$fc))
            if (sum(EtoS[, 2] == i) == 0) {
                resType[seq_len(length(resType))] <- -1
                resMat[seq_len(length(resType))] <- -1
            } else {
                if (sum(rownames(NEMlist$fc) %in%
                        names(which(EtoS[, 2] == i))) == 1) {
                    data.tmp <-
                        t(as.matrix(
                            NEMlist$fc[which(rownames(NEMlist$fc) %in%
                                             names(which(EtoS[, 2] == i))), ]))
                    rownames(data.tmp) <-
                        rownames(NEMlist$fc)[
                            which(rownames(NEMlist$fc) %in%
                                  names(which(EtoS[, 2] == i)))]
                } else {
                    data.tmp <-
                        NEMlist$fc[which(rownames(NEMlist$fc) %in%
                                         names(which(EtoS[, 2] == i))), ]
                }
                resVec <- sum(abs(cor(SCompMat[i, ], t(data.tmp),
                                      method = method)))
                for (j in seq_len(ncol(data.tmp))) { # parallel this!

                    if (verbose) {
                        cat('\r',
                            paste(floor(((i-1)*ncol(data.tmp) + j)/
                                        (ncol(CNOlist@signals[[1]])*
                                         ncol(data.tmp))*100), "%",
                                  sep = ""))
                        flush.console()
                    }

                    sgene <- SCompMat[i, ]
                    mem <- sgene[j]
                    if (mem == 0) {
                        resType[j] <- 1
                        resType[(j+ncol(data.tmp))] <- 1
                        sgene[j] <- 1
                        resMat[j] <- sum(abs(cor(sgene, t(data.tmp),
                                                 method = method)))
                        sgene[j] <- -1
                        resMat[(j+ncol(data.tmp))] <-
                            sum(abs(cor(sgene, t(data.tmp),
                                        method = method)))
                    }
                    if (mem == 1) {
                        resType[j] <- -1
                        resType[(j+ncol(data.tmp))] <- 1
                        sgene[j] <- 0
                        resMat[j] <- sum(abs(cor(sgene, t(data.tmp),
                                                 method = method)))
                        sgene[j] <- -1
                        resMat[(j+ncol(data.tmp))] <-
                            sum(abs(cor(sgene,
                                        t(data.tmp),
                                        method = method)))
                    }
                    if (mem == -1) {
                        resType[j] <- 1
                        resType[(j+ncol(data.tmp))] <- -1
                        sgene[j] <- 1
                        resMat[j] <- sum(abs(cor(sgene, t(data.tmp),
                                                 method = method)))
                        sgene[j] <- 0
                        resMat[(j+ncol(data.tmp))] <-
                            sum(abs(cor(sgene,
                                        t(data.tmp), method = method)))
                    }
                }
            }
            return(list(resMat = resMat, resType = resType, resVec = resVec))
        }

        if (!is.null(parallel)) {
            if (is.list(parallel)) {
                if (length(parallel[[1]]) != length(parallel[[2]])) {
                    stop("The nodes (second list object in parallel) and the
number of cores used on every node (first list object in parallel) must be the
same.") }
                hosts <- character()
                for (i in seq_len(length(parallel[[1]]))) {
                    hosts <- c(hosts, rep(parallel[[2]][i], parallel[[1]][i]))
                }
                hosts <- as.list(hosts)
                sfInit(parallel=TRUE, socketHosts=hosts)
            } else {
                sfInit(parallel=TRUE, cpus=parallel, type = "SOCK")
            }
            sfLibrary("CellNOptR")
        }

        if (!is.null(parallel)) {
            resTmp <- sfLapply(as.list(seq_len(nrow(resMat))), checkSgene)
            sfStop()
        } else {
            resTmp <- lapply(as.list(seq_len(nrow(resMat))), checkSgene)
        }

        for (i in seq_len(nrow(resMat))) {
            resMat[i, ] <- resTmp[[i]]$resMat
            resType[i, ] <- resTmp[[i]]$resType
            resVec[i] <- resTmp[[i]]$resVec[1]
        }
        resType <- resType*(-1)
        resDiff <- (resVec - resMat)/nrow(NEMlist$fc)
        colnames(resDiff) <- c(colnames(NEMlist$fc), colnames(NEMlist$fc))
        rownames(resDiff) <- colnames(CNOlist@signals[[1]])
        resDiff1 <- cbind(resDiff[, seq_len(ncol(NEMlist$fc))], max(resDiff),
                          resDiff[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])

        p1 <- HeatmapOP(resDiff1, Rowv = FALSE, Colv = FALSE, main = main,
                        sub = sub, bordercol = "grey", ...)

        resDiff2 <- cbind(resDiff[, seq_len(ncol(NEMlist$fc))], min(resDiff),
                          resDiff[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])
        resType2 <- cbind(resType[, seq_len(ncol(NEMlist$fc))], min(resType),
                          resType[, (ncol(NEMlist$fc)+1):(2*ncol(NEMlist$fc))])
        resDiff2[which(resDiff2 > 0)] <- 0

        p2 <- HeatmapOP(resDiff2, Colv = FALSE, Rowv = FALSE, main = main,
                        sub = sub, bordercol = "grey", ...)

        resDiff3 <- resDiff2*resType2

        p3 <- HeatmapOP(resDiff3, Colv = FALSE, Rowv = FALSE, main = main,
                        sub = sub, bordercol = "grey", ...)
        res.breaks <-
            seq(-max(abs(min(resDiff3, na.rm = TRUE)),
                     abs(max(resDiff3, na.rm = TRUE))),
                max(abs(min(resDiff3, na.rm = TRUE)), abs(max(resDiff3,
                                                              na.rm = TRUE))),
                (max(abs(min(resDiff3, na.rm = TRUE)), abs(max(resDiff3,
                                                               na.rm = TRUE))) -
                 -max(abs(min(resDiff3, na.rm = TRUE)),
                      abs(max(resDiff3,
                              na.rm = TRUE))))/100)

        p1 <- HeatmapOP(resDiff3[, seq_len(ncol(NEMlist$fc))],
                        bordercol = "grey", Colv = FALSE, Rowv = FALSE,
                        main = "residuals (positive effects)", sub = "",
                        xrot = "60", breaks = res.breaks, colorkey = FALSE)

        p2 <- HeatmapOP(resDiff3[, (ncol(NEMlist$fc)+2):(2*ncol(NEMlist$fc)+1)],
                        bordercol = "grey", Colv = FALSE, Rowv = FALSE,
                        main = "residuals (negative effects)", sub = "",
                        xrot = "60", breaks = res.breaks, colorkey = TRUE)

        if (verbose) {
            print(p1, position=c(0, 0, .48, 1), more=TRUE)
            print(p2, position=c(.48, 0, 1, 1))
        }
        if (cut & all(is.na(resDiff) == FALSE)) {
            if (sum(apply(abs(resDiff1), 1, sum) == 0) > 0) {
                resDiff1 <-
                    resDiff1[-which(apply(abs(resDiff1), 1, sum) == 0), ]
            }
            if (sum(apply(abs(resDiff1), 2, sum) == 0) > 0) {
                resDiff1 <-
                    resDiff1[, -which(apply(abs(resDiff1), 2, sum) == 0)]
            }
            if (sum(apply(abs(resDiff2), 1, sum) == 0) > 0) {
                resDiff2 <-
                    resDiff2[-which(apply(abs(resDiff2), 1, sum) == 0), ]
            }
            if (sum(apply(abs(resDiff2), 2, sum) == 0) > 0) {
                resDiff2 <-
                    resDiff2[, -which(apply(abs(resDiff2), 2, sum) == 0)]
            }
            if (sum(apply(abs(resDiff3), 1, sum) == 0) > 0) {
                resDiff3 <-
                    resDiff3[-which(apply(abs(resDiff3), 1, sum) == 0), ]
            }
            if (sum(apply(abs(resDiff3), 2, sum) == 0) > 0) {
                resDiff3 <-
                    resDiff3[, -which(apply(abs(resDiff3), 2, sum) == 0)]
            }
        }
        return(list(resDiff1 = resDiff1, resDiff2 = resDiff2,
                    resDiff3 = resDiff3))
    }
#' reduce graph
#'
#' reduces the size of a graph if possible to an equivalent sub-graph
#' @param bString binary vector indicating the sub-graph given a model
#' @param model model object for the whole graph space
#' @param CNOlist CNOlist object
#' @author Martin Pirkl
#' @return equivalent sub-graph denoted by a bString
#' @export
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signal = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' bString <- reduceGraph(rep(1, length(model$reacID)), model, CNOlist)
reduceGraph <-
    function(bString, model, CNOlist) {
        if (any(bString != 0)) {
            stimuli <- colnames(CNOlist@stimuli)
            graph <- model$reacID[which(bString == 1)]
            tmp <- unlist(strsplit(graph, "="))
            tmp <- unlist(strsplit(tmp, "\\+"))
            tmp <- unique(gsub("!", "", tmp))
            for (i in tmp) {
                if (!(i %in% stimuli) &
                    length(grep(paste("=", i, sep = ""), graph)) == 0) {
                    get <- grep(paste("^", i, sep = ""), graph)
                    if (length(get) > 0) {
                        graph <- graph[-get]
                    }
                }
            }
            bString <- numeric(length(bString))
            bString[which(model$reacID %in% graph)] <- 1
        }
        bString <- absorption(bString, model)
        return(bString)
    }
#' simulate states
#'
#' simulates the activation pattern (truth table) of a hyper-graph and
#' annotated perturbation experiments
#' @param CNOlist CNOlist object
#' @param model model object
#' @param bString binary vector denoting the sub-graph given model
#' @param NEMlist NEMlist object only for devel
#' @author Martin Pirkl
#' @return return the truth tables for certain perturbation experiments
#' as a numeric matrix
#' @export
#' @import
#' CellNOptR
#' matrixStats
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signal = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' states <- simulateStatesRecursive(CNOlist, model,
#' rep(1, length(model$reacID)))
simulateStatesRecursive <-
    function(CNOlist, model, bString, NEMlist = NULL) {
        getState <- function(CNOlist, node, signalStates, graph,
                             children = NULL, NEMlist = NULL) {
            graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
            if (length(graphCut) == 0) {
                signalStates[, node] <- 0
            } else {
                sop <- numeric(nrow(signalStates))
                children2 <- gsub("!", "", children)
                for (i in graphCut) {
                    parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
                    if (length(parents) == 0) {
                        pob <- rep(0, nrow(signalStates))
                    } else {
                        pob <- rep(1, nrow(signalStates))
                    }
                    for (j in parents) {
                        j2 <- gsub("!", "", j)
                        if (any(is.na(signalStates[, j2]) == TRUE)) {
                            if (j %in% j2) {
                                node2 <- node
                                add1 <- 0
                            } else {
                                node2 <- paste("!", node, sep = "")
                                add1 <- 1
                            }
                            if (j2 %in% children2) {
                                ## this speeds up the process and will
                                subGraph <- graph
                                subGraph <-
                                    subGraph[
                                        -grep(paste(".*=", node, "|.*", j2,
                                                    ".*=.*", sep = ""),
                                              subGraph)]
                                signalStatesTmp <-
                                    getState(CNOlist = CNOlist,
                                             node = j2,
                                             signalStates = signalStates,
                                             graph = subGraph,
                                             children = children2[
                                                 -which(children2 %in% node)],
                                             NEMlist)
                                if (add1 == 0) {
                                    pobMult <- signalStatesTmp[, j2]
                                } else {
                                    pobMult <- add1 - signalStatesTmp[, j2]
                                }
                            } else {
                                signalStates <-
                                    getState(CNOlist = CNOlist, node = j2,
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
                if (node %in% colnames(CNOlist@inhibitors)) {
                    sop <- sop*(1 - CNOlist@inhibitors[, node])
                }
                if (node %in% colnames(CNOlist@stimuli)) {
                    sop <- max(sop, CNOlist@stimuli[, node])
                }
                signalStates[, node] <- sop
            }
            return(signalStates)
        }
        bString <- reduceGraph(bString, model, CNOlist)
        stimuli <- colnames(CNOlist@stimuli)
        signals <-
            sort(c(colnames(CNOlist@inhibitors),
                   model$namesSpecies[
                             -which(model$namesSpecies %in%
                                    c(stimuli,
                                      colnames(CNOlist@inhibitors)))]))
        graph0 <- model$reacID[which(bString == 1)]
        stimuliStates <- CNOlist@stimuli
        if (!is.null(NEMlist$signalStates)) {
            signalStates <- NEMlist$signalStates
        } else {
            signalStates <- matrix(NA, nrow = nrow(CNOlist@signals[[2]]),
                                   ncol = length(signals))
            rownames(signalStates) <- rownames(CNOlist@signals[[2]])
            colnames(signalStates) <- signals
            signalStates <- cbind(stimuliStates, signalStates)
        }
        for (k in signals) {
            if (sum(is.na(signalStates[, k]) == TRUE) ==
                length(signalStates[, k])) {
                signalStates <- getState(CNOlist = CNOlist, node = k,
                                         signalStates = signalStates,
                                         graph = graph0, children = NULL,
                                         NEMlist)
            }
        }
        signalStates <- signalStates[, which(colnames(signalStates) %in%
                                             colnames(CNOlist@signals[[1]]))]
        if (ncol(CNOlist@signals[[1]]) != 1) {
            signalStates <- signalStates[, order(colnames(signalStates))]
        } else {
            signalStates <- as.matrix(signalStates)
            colnames(signalStates) <- colnames(CNOlist@signals[[1]])
        }
        return(signalStates = signalStates)
    }
#' transitive closure
#'
#' calculates transitive closure of a hyper-graph
#' @param g hyper-graph in normal form
#' @param max.iter maximal iteration till convergence
#' @param verbose verbose output?
#' @author Martin Pirkl
#' @return transitive closure in normal form
#' @export
#' @examples
#' g <- c("A=B", "B=C")
#' gclose <- transClose(g)
transClose <-
    function(g, max.iter = NULL, verbose = FALSE) {
        v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")),
                                                  "\\+"))))
        if (is.null(max.iter)) {
            h <- getHierarchy(g)
            max.iter <- length(h) - 2 # should be sufficient
        }
        if (verbose) {
            print(paste("maximum iterations: ", max.iter, sep = ""))
        }
        g.out <- unique(gsub(".*=", "", g))
        g.closed <- g
        for (iter in seq_len(max.iter)) {
            g.old <- g.closed

            if (verbose) {
                cat('\r', paste("iteration: ", iter, sep = ""))
                flush.console()
            }
            for (i in g.closed) {
                input <-
                    gsub("!", "", unlist(strsplit(unlist(strsplit(i, "="))[1],
                                                  "\\+")))
                input <- intersect(input, g.out)
                output <- gsub(".*=", "", i)
                if (length(input) == 0) { next() }
                for (j in unique(input)) {
                    if (j %in% unlist(strsplit(unlist(strsplit(i, "="))[1],
                                               "\\+"))) {
                        for (k in gsub("=.*", "",
                                       g[grep(paste("=", j, sep = ""),
                                              g)])) {
                            g.closed <- c(g.closed, gsub(j, k, i))
                        }
                    } else {
                        literals <- list()
                        count <- 0
                        for (k in unique(gsub("=.*", "",
                                              g[grep(paste("=", j,
                                                           sep = ""), g)]))) {
                            count <- count + 1
                            literals[[count]] <-
                                gsub("!!", "", paste("!",
                                                     unlist(strsplit(k, "\\+")),
                                                     sep = ""))
                        }
                        combis <- expand.grid(literals)
                        combis <- apply(combis, c(1,2), as.character)
                        for (k in seq_len(nrow(combis))) {
                            g.closed <-
                                c(g.closed, gsub(paste("!", j, sep = ""),
                                                 paste(combis[k, ],
                                                       collapse = "+"), i))
                        }
                    }
                }
            }
            g.closed <- unique(g.closed)
            if (all(g.closed %in% g.old)) {
                if (verbose) {
                    cat("\n")
                    print(paste("successfull convergence", sep = ""))
                }
                break()
            }
        }
        if (verbose) {
            cat("\n")
        }
        g.closed <- unique(g.closed)
        for (j in seq_len(length(g.closed))) {
            i <- g.closed[j]
            input <- unlist(strsplit(i, "="))
            output <- input[2]
            input <- unlist(strsplit(input[1], "\\+"))
            input <- unique(input)
            g.closed[j] <- paste(paste(input, collapse = "+"),
                                 output, sep = "=")
        }
        if (length(grep(paste(paste(v, ".*", v, ".*=", sep = ""),
                              collapse = "|"), g.closed)) > 0) {
            g.closed <- g.closed[-grep(paste(paste(v, ".*", v, ".*=", sep = ""),
                                             collapse = "|"), g.closed)]
        }
        return(g.closed)
    }
#' transitive reduction
#'
#' calculates transitive reduction of a hyper-graph in normal form
#' @param g hyper-graph in normal form
#' @param max.iter maximal number of iterations till convergence
#' @param verbose verbose output?
#' @author Martin Pirkl
#' @return transitive reduction of the hyper-graph in normal form
#' @export
#' @examples
#' g <- c("A=B", "A=C", "B=C", "B=D", "!A=D")
#' gred <- transRed(g)
transRed <-
    function(g, max.iter = NULL, verbose = FALSE) {
        v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")),
                                                  "\\+"))))
        if (is.null(max.iter)) {
            max.iter <- length(v) - 2
        }
        a <- dnf2adj(g)
        g2 <- g
        h <- getHierarchy(g2)
        if (length(h) > 2) {
            for (i in seq_len((length(h)-2))) {
                for (j in h[[i]]) {
                    for (k in (i+2):length(h)) {
                        for (l in h[[k]]) {
                            if (length(grep(paste(".*", j,
                                                  ".*=", l, sep = ""), g2)) !=
                                0) {
                                if (length(grep(paste(".*=", l,
                                                      sep = ""), g2)) >
                                    1) {
                                    g2 <- g2[-grep(paste(".*", j, ".*=", l,
                                                         sep = ""), g2)]
                                }
                            }
                        }
                    }
                }
            }
        }
        g3 <- transClose(g2, max.iter)
        if (sum(g %in% g3) > 0) {
            g4 <- g[-which(g %in% g3)]
        }
        g5 <- unique(c(g2, g4))
        return(g5)
    }
#' validate graph
#'
#' plotting the observed response scheme of an effect reporter and the
#' expected response scheme of the regulating signalling gene
#' @param CNOlist CNOlist object.
#' @param fc matrix of foldchanges (observed response scheme, ORS).
#' @param exprs optional matrix of normalized expression (observed
#' activation scheme).
#' @param approach default "fc" for foldchanges or signed effect probabilities,
#' "abs" for absolute effects or probabilities
#' @param model Model object.
#' @param bString Binary string denoting the hyper-graph.
#' @param Egenes Maximal number of visualized E-genes.
#' @param Sgene Integer denoting the S-gene. See
#' colnames(CNOlist@signals[[1]]) to match integer with S-gene name.
#' @param parameters see ?bnem
#' @param plot Plot the heatmap. If FALSE, only corresponding
#' information is printed.
#' @param disc Discretize the data.
#' @param affyIds Experimental. Turn Affymetrix Ids into HGNC
#' gene symbols.
#' @param relFit see ?bnem
#' @param xrot See package ?epiNEM::HeatmapOP
#' @param Rowv See package ?epiNEM::HeatmapOP
#' @param Colv See package ?epiNEM::HeatmapOP
#' @param dendrogram See package ?epiNEM::HeatmapOP
#' @param soft weighting the expected pattern, if TRUE
#' @param colSideColors See package ?epiNEM::HeatmapOP
#' @param affychip Define Affymetrix chip used to generate the data
#' (optional).
#' @param method see ?bnem
#' @param ranks Turns data into ranks, if TRUE
#' @param breaks See package ?epiNEM::HeatmapOP
#' @param col See package ?epiNEM::HeatmapOP
#' @param sizeFac size factor penelizing the hyper-graph size.
#' @param verbose verbose output, if TRUE
#' @param order Order by "rank" or "name" or "none"
#' @param ... additional arguments for ?epiNEM::HeatmapOP
#' @author Martin Pirkl
#' @return lattice object with matrix information
#' @export
#' @import
#' CellNOptR
#' stats
#' RColorBrewer
#' epiNEM
#' matrixStats
#' @examples
#' sifMatrix <- rbind(c("A", 1, "B"), c("A", 1, "C"), c("B", 1, "D"),
#' c("C", 1, "D"))
#' write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
#' col.names = FALSE,
#' quote = FALSE)
#' PKN <- CellNOptR::readSIF("temp.sif")
#' unlink('temp.sif')
#' CNOlist <- dummyCNOlist("A", c("B","C","D"), maxStim = 1, maxInhibit = 2,
#' signal = c("A", "B","C","D"))
#' model <- CellNOptR::preprocessing(CNOlist, PKN, maxInputsPerGate = 100)
#' exprs <- matrix(rnorm(nrow(slot(CNOlist, "cues"))*10), 10,
#' nrow(slot(CNOlist, "cues")))
#' fc <- computeFc(CNOlist, exprs)
#' initBstring <- rep(0, length(model$reacID))
#' res <- bnem(search = "greedy", CNOlist = CNOlist, fc = fc,
#' model = model, parallel = NULL, initBstring = initBstring, draw = FALSE,
#' verbose = FALSE, maxSteps = Inf)
#' rownames(fc) <- seq_len(nrow(fc))
#' val <- validateGraph(CNOlist = CNOlist, fc = fc, model = model,
#' bString = res$bString, Egenes = 10, Sgene = 4)
validateGraph <-
    function(CNOlist, fc=NULL, exprs=NULL, approach = "fc", model, bString,
             Egenes = 25, Sgene = 1,
             parameters = list(cutOffs = c(0,1,0), scoring = c(0.1,0.2,0.9)),
             plot = TRUE,
             disc = 0, affyIds = TRUE, relFit = FALSE,
             xrot = 25, Rowv = FALSE, Colv = FALSE,
             dendrogram = "none", soft = TRUE, colSideColors = NULL,
             affychip = "hgu133plus2", method = "s", ranks = FALSE,
             breaks = NULL, col = "RdYlGn", sizeFac = 10^-10,
             verbose = TRUE, order = "rank", ...) {
        csc <- TRUE
        colnames <- "bio"
        complete <- FALSE
        sim <- 0

        NEMlist <- list()
        NEMlist$fc <- fc
        NEMlist$exprs <- exprs

        myCN2bioCN <- function(x, stimuli, inhibitors) {
            y <- gsub("_vs_", ") vs (", x)
            for (i in inhibitors) {
                y <- gsub(i, paste(i, "\\-", sep = ""), y)
            }
            for (i in stimuli) {
                y <- gsub(i, paste(i, "\\+", sep = ""), y)
            }
            y <- gsub("Ctrl", "control", paste("(", gsub("_", ",", y), ")",
                                               sep = ""))
            return(y)
        }

        gene2protein <- function(genes, strict = FALSE) {
            if (strict) {
                gene2prot <- cbind(
                    c("^APC$", "^ATF2$", "^BIRC2$", "^BIRC3$", "^CASP4$",
                      "^CASP8$", "^CFLAR$", "^CHUK$", "^CTNNB1$", "^DKK1$",
                      "^DKK4$", "^FLASH$", "^IKBKB$", "^IKBKG$", "^JUN$",
                      "^MAP2K1$", "^MAP3K14$", "^MAP3K7$", "^MAPK8$",
                      "^PIK3CA$", "^RBCK1$", "^RELA$", "^RIPK1$", "^RIPK3$",
                      "^RNF31$", "^SHARPIN$", "^TAB2$", "^TCF4$", "^TCF7L2$",
                      "^TNFRSF10A$", "^TNFRSF10B$", "^TNFRSF1A$", "^TNIK$",
                      "^TRAF2$", "^USP2$", "^WLS$", "^WNT11$", "^WNT5A$",
                      "^TNFa$", "^TRAIL"),
                    c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8",
                      "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4",
                      "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik",
                      "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1",
                      "Rip3", "Hoip", "Sharpin", "Tab2", "fake", "Tcf4",
                      "Dr4", "Dr5", "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi",
                      "Wnt11", "Wnt5A", "Tnfa", "Trail")
                )
            } else {
                gene2prot <- cbind(
                    c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8",
                      "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH",
                      "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14", "MAP3K7",
                      "MAPK8", "PIK3CA", "RBCK1", "RELA", "RIPK1", "RIPK3",
                      "RNF31", "SHARPIN", "TAB2", "TCF4", "TCF7L2",
                      "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "TNIK", "TRAF2",
                      "USP2", "WLS", "WNT11", "WNT5A", "TNFa", "TRAIL"),
                    c("Apc", "Atf2", "cIap1", "cIap2", "Casp4", "Casp8",
                      "c-Flip", "Ikka", "Beta-Cat.", "Dkk1", "Dkk4",
                      "Casp8ap2", "Ikkb", "Nemo", "cJun", "Mekk", "Nik",
                      "Tak1", "Jnk", "Pi3k", "Hoil1", "RelA", "Rip1", "Rip3",
                      "Hoip", "Sharpin", "Tab2", "fake", "Tcf4", "Dr4", "Dr5",
                      "Tnfr1", "Tnik", "Traf2", "Usp2", "Evi", "Wnt11",
                      "Wnt5A", "Tnfa", "Trail")
                )
            }
            for (i in seq_len(nrow(gene2prot))) {
                genes <- gsub(gene2prot[i, 1], gene2prot[i, 2], genes)
            }
            return(genes)
        }

        protein2gene <- function(proteins, strict = FALSE) {
            if (strict) {
                gene2prot <- cbind(
                    c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8",
                      "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH",
                      "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14",
                      "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA",
                      "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2",
                      "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B",
                      "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS",
                      "WNT11", "WNT5A", "TNFa", "TRAIL"),
                    c("^Apc$", "^Atf2$", "^cIap1$", "^cIap2$",
                      "^Casp4$", "^Casp8$", "^c-Flip$", "^Ikka$",
                      "^Beta-Cat.$", "^Dkk1$", "^Dkk4$", "^Casp8ap2$",
                      "^Ikkb$", "^Nemo$", "^cJun$", "^Mekk$", "^Nik$",
                      "^Tak1$", "^Jnk$", "^Pi3k$", "^Hoil1$", "^RelA$",
                      "^Rip1$", "^Rip3$", "^Hoip$", "^Sharpin$", "^Tab2$",
                      "^fake$", "^Tcf4$", "^Dr4$", "^Dr5$", "^Tnfr1$",
                      "^Tnik$", "^Traf2$", "^Usp2$", "^Evi$", "^Wnt11$",
                      "^Wnt5A$", "^Tnfa$", "^Trail$")
                )
            } else {
                gene2prot <- cbind(
                    c("APC", "ATF2", "BIRC2", "BIRC3", "CASP4", "CASP8",
                      "CFLAR", "CHUK", "CTNNB1", "DKK1", "DKK4", "FLASH",
                      "IKBKB", "IKBKG", "JUN", "MAP2K1", "MAP3K14",
                      "MAP3K7", "MAPK8", "PIK3CA", "RBCK1", "RELA",
                      "RIPK1", "RIPK3", "RNF31", "SHARPIN", "TAB2",
                      "TCF4", "TCF7L2", "TNFRSF10A", "TNFRSF10B",
                      "TNFRSF1A", "TNIK", "TRAF2", "USP2", "WLS",
                      "WNT11", "WNT5A", "TNFa", "TRAIL"),
                    c("Apc", "Atf2", "cIap1", "cIap2", "Casp4",
                      "Casp8", "c-Flip", "Ikka", "Beta-Cat.", "Dkk1",
                      "Dkk4", "Casp8ap2", "Ikkb", "Nemo", "cJun",
                      "Mekk", "Nik", "Tak1", "Jnk", "Pi3k", "Hoil1",
                      "RelA", "Rip1", "Rip3", "Hoip", "Sharpin", "Tab2",
                      "fake", "Tcf4", "Dr4", "Dr5", "Tnfr1", "Tnik",
                      "Traf2", "Usp2", "Evi", "Wnt11", "Wnt5A", "Tnfa",
                      "Trail")
                )
            }
            for (i in seq_len(nrow(gene2prot))) {
                proteins <- gsub(gene2prot[i, 2], gene2prot[i, 1], proteins)
            }
            return(proteins)
        }

        colSideColorsSave <- NULL
        bad.data <- FALSE
        errorMat <- function() {
            error.mat <- matrix(0, nrow = 50, ncol = (7*3+8))
            error.mat[seq_len(10), c(2:4, 6:8, 10:11,
                                     14:15, 18:20, 22:24, 26:28)] <- 3.5
            error.mat[11:20, c(2,4, 6,8, 10,12,
                               14,16, 18,20, 23, 26,28)] <- 3.5
            error.mat[21:30, c(2:4, 6:8, 10,12, 14,16,
                               18:20, 23, 26:28)] <- 3.5
            error.mat[31:40, c(2,4, 6,8, 10,12, 14,16,
                               18,20, 23, 26,28)] <- 3.5
            error.mat[41:50, c(2:4, 6,8, 10:11, 14:15,
                               18,20, 23, 26,28)] <- 3.5
            error.mat <- cbind(error.mat[, seq_len(13)],
                               matrix(0, nrow = 50, ncol = 5),
                               error.mat[, 14:29])
            error.mat <- rbind(matrix(0, nrow = 10,
                                      ncol = ncol(error.mat)),
                               error.mat, matrix(0, nrow = 10,
                                                 ncol = ncol(error.mat)),
                               error.mat, matrix(0, nrow = 10,
                                                 ncol = ncol(error.mat)),
                               error.mat, matrix(0, nrow = 10,
                                                 ncol = ncol(error.mat)),
                               error.mat, matrix(0, nrow = 10,
                                                 ncol = ncol(error.mat)))
            rownames(error.mat) <- seq_len(250)
            colnames(error.mat) <- seq_len(((7*3+8)+5))
            error.mat <- error.mat
            error.mat <- error.mat
            error.mat <- error.mat[rep(seq_len(nrow(error.mat)), each = 1),
                                   rep(seq_len(ncol(error.mat)), each = 5)]
            error.mat <- smoothMatrix(error.mat, 5, torus = FALSE)
            return(error.mat)
        }
        CNOlist <- checkCNOlist(CNOlist)
        NEMlist <- checkNEMlist(NEMlist = NEMlist, CNOlist = CNOlist,
                                parameters = parameters, approach = approach,
                                method = method)
        tmp <- computeScoreNemT1(CNOlist, model = model, bString, NAFac = 1,
                                 approach = approach, NEMlist = NEMlist,
                                 tellme = 1, parameters = parameters,
                                 sim = sim, relFit = relFit, method = method,
                                 sizeFac = sizeFac, verbose = FALSE)
        EtoS <- tmp$EtoS
        if (verbose) {
            print(paste(Sgene, ".", colnames(CNOlist@signals[[1]])[Sgene],
                        ": ", sum(EtoS[, 2] == Sgene), sep = ""))
            print(paste("Activated: ", sum(EtoS[, 2] == Sgene & EtoS[, 3] == 1),
                        sep = ""))
            print(paste("Inhibited: ", sum(EtoS[, 2] == Sgene &
                                           EtoS[, 3] == -1), sep = ""))
            print("Summary Score:")
            print(summary(EtoS[which(EtoS[, 2] == Sgene), 4]))
            dups <- sum(duplicated(rownames(EtoS)) == TRUE)
            if (dups > 0) {
                used <- length(unique(rownames(EtoS)))
            } else {
                used <- nrow(EtoS)
            }
            print(paste("Unique genes used: ", (used), " (",
                        round((used/nrow(NEMlist$fc))*100, 2), " %)",
                        sep = ""))
            print(paste("Duplicated genes: ", dups, sep = ""))
            print("Overall fit:")
            print(summary(EtoS[, 4]))
        }
        indexList <- NULL
        if (is.null(indexList) == TRUE) {
            indexList = indexFinder(CNOlist, model)
        }
        modelCut = cutModel(model, bString)
        if (sim == 0) {
            simResults <-
                simulateStatesRecursive(CNOlist = CNOlist,
                                        model = modelCut,
                                        bString =
                                            (numeric(length(modelCut$reacID)) +
                                             1))
        }
        if (sim == 1) {
            simResults <-
                simulateStatesRecursiveAdd(CNOlist = CNOlist, model = modelCut,
                                           bString =
                                               (numeric(
                                               length(modelCut$reacID)) + 1))
        }
        rownames(simResults) <- seq_len(nrow(simResults))
        simResults <- simResults[, which(colnames(simResults) %in%
                                         colnames(CNOlist@signals[[1]]))]
        SCompMat <- computeFc(CNOlist, t(simResults))
        SCompMat <- SCompMat[, colnames(NEMlist$fc)]

        if (parameters$cutOffs[3] == -1) {
            method <- checkMethod(method)
            S.mat <- SCompMat
            ## do median polish over gene clusters
            data.med <- NEMlist$fc[seq_len(ncol(CNOlist@signals[[2]])), ]*0
            Epos <- EtoS[order(rownames(EtoS)), seq_len(2)]
            for (i in seq_len(ncol(CNOlist@signals[[1]]))) {
                tmp <-
                    medpolish(rbind(
                        NEMlist$fc[Epos[which(Epos[, 2] == i),
                                        1], ],
                        NEMlist$fc[Epos[which(Epos[, 2] ==
                                              (i+ncol(CNOlist@signals[[2]]))),
                                        1], ]), trace.iter=FALSE)
                data.med[i, ] <- tmp$col
            }
            E.mat <- data.med
            E.mat[is.na(E.mat)] <- 0
            tmp <- which(apply(E.mat, 1, sum) != 0)
            E.mat <- E.mat[which(apply(E.mat, 1, sum) != 0), ]
            rownames(E.mat) <- rownames(S.mat)[tmp]
            NEMlist$fc <- E.mat
            if (parameters$scoring[1] == 0) {
                S.mat[which(S.mat == 0)] <- NA
            }
            if ("pearson" %in% method) {
                cosine.sim <- -t(cor(t(S.mat), t(E.mat), method = "p",
                                     use = "pairwise.complete.obs"))
            }
            if ("spearman" %in% method) {
                S.mat.ranks <- apply(S.mat, 1, rank)
                cosine.sim <- -t(cor(S.mat.ranks, t(E.mat), method = "p",
                                     use = "pairwise.complete.obs"))
            }
            cosine.sim[is.na(cosine.sim)] <- 0
            MSEAfc <- cosine.sim
            MSEIfc <- -cosine.sim
            R <- cbind(MSEAfc, MSEIfc)
            R[is.na(R)] <- max(R[!is.na(R)])
            MSEE <- matrixStats::rowMins(R)
        }

        ## matrix visualisation for egenes fitted:
        if ("fc" %in% approach) {
            check.data <- NEMlist$fc
            check.model <- SCompMat
            rownames(check.model) <- colnames(simResults)
        }
        if ("abs" %in% approach) {
            check.data <- NEMlist$exprs # here the $norm is needed
            check.model <- simResults
            colnames(check.model) <- colnames(CNOlist@signals[[2]])
        }

        Egenes <- Egenes

        Egenes <- min(Egenes, sum(EtoS[, 2] == Sgene))

        if (Egenes == 0) {
            mainlab <- paste("Regulated by ",
                             rownames(check.model)[Sgene], "\n", sep = "")
            if (plot) {
                print(HeatmapOP(errorMat(), main = mainlab, sub = "",
                                Colv = FALSE, Rowv = FALSE, col = "RdBu",
                                coln = 12, breaks = seq(0,8, 0.1)), ...)
            }
            genesInfo <- as.matrix(0)
            rownames(genesInfo) <- "dummy"
            return(list(genesInfo = genesInfo, data = genesInfo))
        }

        if (complete) {
            egenefit <- matrix(0, nrow = (sum(EtoS[, 2] == Sgene)+1),
                               ncol = ncol(check.data))
        } else {
            egenefit <- matrix(0, nrow = (Egenes+1), ncol = ncol(check.data))
        }

        egenefit[1,] <- check.model[Sgene, ]
        rownames(egenefit) <- seq_len(nrow(egenefit))
        rownames(egenefit)[1] <- rownames(check.model)[Sgene]
        colnames(egenefit) <- seq_len(ncol(egenefit))
        colnames(egenefit) <- colnames(check.data)
        if (complete) {
            activatedEgenes <- numeric(sum(EtoS[, 2] == Sgene)+1)
        } else {
            activatedEgenes <- numeric(Egenes+1)
        }

        count <- 0
        for (i in seq_len(nrow(EtoS))) {
            if (EtoS[i, 2] == Sgene) {
                egenefit[count+2,] <-
                    check.data[which(rownames(check.data) %in%
                                     rownames(EtoS)[i]), ]
                rownames(egenefit)[count+2] <- rownames(EtoS)[i]
                if (EtoS[i, 3] == 1) {
                    activatedEgenes[count+2] <- 1
                }
                count <- count + 1
            }
            if (count >= Egenes & !complete) { break() }
        }

        Egenes <- count

        if (affyIds == FALSE) {
            temp <-
                as.vector(unlist(mget(unique(rownames(egenefit)[-1]),
                                      get(paste(affychip, "SYMBOL",
                                                sep = "")))))
            temp2 <- rownames(egenefit)[-1]
            if (sum(is.na(temp) == TRUE) > 0) {
                temp[is.na(temp)] <- temp2[is.na(temp)]
            }
            rownames(egenefit)[-1] <- temp
        }

        rownames(egenefit)[is.na(rownames(egenefit))] <- "NA"
        count <- 0
        if (min(egenefit) != max(egenefit)) {
            if (disc == 1) {
                egenefit[1, ] <- disc(egenefit[1, ], 0.5)
                for (i in 2:nrow(egenefit)) {
                    if (parameters$cutOffs[2] == 0) {
                        save.model <- egenefit[1, ]
                        egenefit[1, ] <- save.model
                        real.breaks <-
                            seq(-max(abs(min(egenefit)),
                                     abs(max(egenefit))),
                                max(abs(min(egenefit)),abs(max(egenefit))),0.1)
                        if (length(real.breaks) > 101) {
                            real.breaks <-
                                real.breaks[
                                    c((floor(length(
                                          real.breaks)/2)-50):
                                      floor(length(real.breaks)/2),
                                    (floor(length(real.breaks)/2)+1):
                                    (floor(length(real.breaks)/2)+50))]
                        }
                    } else {
                        if (activatedEgenes[i] == 1) {
                            onePosI <-
                                c(which(egenefit[1, ] == 1 & egenefit[i, ] >
                                        parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 0 & egenefit[i, ] >
                                        parameters$cutOffs[2]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] >
                                        parameters$cutOffs[2]))
                            onePosII <-
                                which(egenefit[1, ] == 1 & egenefit[i, ] >
                                      parameters$cutOffs[1] & egenefit[i, ] <=
                                      parameters$cutOffs[2])
                            oneNegI <-
                                c(which(egenefit[1, ] == -1 & egenefit[i, ] <
                                        -parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 0 & egenefit[i, ] <
                                        -parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] <
                                        -parameters$cutOffs[2]))
                            oneNegII <-
                                which(egenefit[1, ] == -1 & egenefit[i, ] <
                                      -parameters$cutOffs[1] & egenefit[i, ] >=
                                      -parameters$cutOffs[2])
                            zeros <-
                                c(which(egenefit[1, ] == 0 &
                                        abs(egenefit[i, ]) <=
                                        parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 1 &
                                        egenefit[i, ] <= parameters$cutOffs[1] &
                                        egenefit[i, ] >=
                                        -parameters$cutOffs[2]),
                                  which(egenefit[1, ] == -1 &
                                        egenefit[i, ] >=
                                        -parameters$cutOffs[1] &
                                        egenefit[i, ] <= parameters$cutOffs[2]))
                            zerosI <-
                                c(which(egenefit[1, ] == 0 &
                                        abs(egenefit[i, ]) <=
                                        parameters$cutOffs[1]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] <=
                                        parameters$cutOffs[1] & egenefit[i, ] >=
                                        -parameters$cutOffs[1]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] <=
                                        parameters$cutOffs[1] & egenefit[i, ] >=
                                        -parameters$cutOffs[1]))
                            zerosII <-
                                c(which(egenefit[1, ] == 0 & egenefit[i, ] <=
                                        parameters$cutOffs[2] & egenefit[i, ] >
                                        parameters$cutOffs[1]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] >
                                        parameters$cutOffs[1] & egenefit[i, ] <=
                                        parameters$cutOffs[2]))
                            zerosIII <-
                                c(which(egenefit[1, ] == 0 & egenefit[i, ] >=
                                        -parameters$cutOffs[2] & egenefit[i, ] <
                                        -parameters$cutOffs[1]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] <
                                        -parameters$cutOffs[1] &
                                        egenefit[i, ] >=
                                        -parameters$cutOffs[2]))
                            if (soft) {
                                egenefit[i, onePosI] <- 3
                                egenefit[i, oneNegI] <- -3
                                egenefit[i, onePosII] <- 2
                                egenefit[i, oneNegII] <- -2
                                egenefit[i, zerosI] <- 0
                                egenefit[i, zerosII] <- 0.5
                                egenefit[i, zerosIII] <- -0.5
                            } else {
                                egenefit[i, onePosI] <- 2
                                egenefit[i, oneNegI] <- -2
                                egenefit[i, onePosII] <- 2
                                egenefit[i, oneNegII] <- -2
                                egenefit[i, zerosI] <- 0
                                egenefit[i, zerosII] <- 0
                                egenefit[i, zerosIII] <- 0
                            }
                        } else {
                            onePosI <-
                                c(which(egenefit[1, ] == 1 & egenefit[i, ] >
                                        parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 0 & egenefit[i, ] >
                                        parameters$cutOffs[2]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] >
                                        parameters$cutOffs[2]))
                            onePosII <-
                                which(egenefit[1, ] == -1 & egenefit[i, ] >
                                      parameters$cutOffs[1] & egenefit[i, ] <=
                                      parameters$cutOffs[2])
                            oneNegI <-
                                c(which(egenefit[1, ] == -1 & egenefit[i, ] <
                                        -parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 0 & egenefit[i, ] <
                                        -parameters$cutOffs[2]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] <
                                        -parameters$cutOffs[2]))
                            oneNegII <-
                                which(egenefit[1, ] == 1 & egenefit[i, ] <
                                      -parameters$cutOffs[1] & egenefit[i, ] >=
                                      -parameters$cutOffs[2])
                            zerosI <-
                                c(which(egenefit[1, ] == 0 &
                                        abs(egenefit[i, ]) <=
                                        parameters$cutOffs[1]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] <=
                                        parameters$cutOffs[1] & egenefit[i, ] >=
                                        -parameters$cutOffs[1]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] <=
                                        parameters$cutOffs[1] & egenefit[i, ] >=
                                        -parameters$cutOffs[1]))
                            zerosII <-
                                c(which(egenefit[1, ] == 0 & egenefit[i, ] <=
                                        parameters$cutOffs[2] & egenefit[i, ] >
                                        parameters$cutOffs[1]),
                                  which(egenefit[1, ] == 1 & egenefit[i, ] >
                                        parameters$cutOffs[1] & egenefit[i, ] <=
                                        parameters$cutOffs[2]))
                            zerosIII <-
                                c(which(egenefit[1, ] == 0 & egenefit[i, ] >=
                                        -parameters$cutOffs[2] & egenefit[i, ] <
                                        -parameters$cutOffs[1]),
                                  which(egenefit[1, ] == -1 & egenefit[i, ] <
                                        -parameters$cutOffs[1] &
                                        egenefit[i, ] >=
                                        -parameters$cutOffs[2]))
                            if (soft) {
                                egenefit[i, onePosI] <- 3
                                egenefit[i, oneNegI] <- -3
                                egenefit[i, onePosII] <- 2
                                egenefit[i, oneNegII] <- -2
                                egenefit[i, zerosI] <- 0
                                egenefit[i, zerosII] <- 0.5
                                egenefit[i, zerosIII] <- -0.5
                            } else {
                                egenefit[i, onePosI] <- 2
                                egenefit[i, oneNegI] <- -2
                                egenefit[i, onePosII] <- 2
                                egenefit[i, oneNegII] <- -2
                                egenefit[i, zerosI] <- 0
                                egenefit[i, zerosII] <- 0
                                egenefit[i, zerosIII] <- 0
                            }
                        }
                        real.breaks <- c(-4.5,-3.5,-2.5,-1.5,-0.75,
                                         -0.25,0.25,0.75,1.5,2.5,3.5,4.5)
                    }
                }
                egenefit[1, ] <- egenefit[1, ]*4
                if (!is.null(colSideColors)) {
                    colSideColorsSave <- colSideColors
                }
                if (min(egenefit) != max(egenefit)) {
                    if (csc) {
                        colSideColors <- rep("black", ncol(egenefit))
                        for (i in seq_len(length(colnames(egenefit)))) {
                            names <-
                                unlist(strsplit(colnames(egenefit)[i], "_"))
                            if (length(names) > 1) {
                                if (names[1] %in%
                                    colnames(CNOlist@inhibitors) &
                                    names[length(names)] %in%
                                    colnames(CNOlist@stimuli)) {
                                    colSideColors[i] <- "brown"
                                }
                                if (names[1] %in%
                                    colnames(CNOlist@stimuli) &
                                    names[length(names)] %in%
                                    colnames(CNOlist@inhibitors)) {
                                    colSideColors[i] <- "orange"
                                }
                                if (names[1] %in% "Ctrl" &
                                    names[length(names)] %in%
                                    colnames(CNOlist@stimuli)) {
                                    colSideColors[i] <- "yellow"
                                }
                                if (names[1] %in% "Ctrl" &
                                    names[length(names)] %in%
                                    colnames(CNOlist@inhibitors)) {
                                    colSideColors[i] <- "blue"
                                }
                            } else {
                                if (names %in% colnames(CNOlist@inhibitors)) {
                                    colSideColors[i] <- "blue"
                                } else {
                                    colSideColors[i] <- "yellow"
                                }
                            }
                        }
                        if (!is.null(colSideColorsSave)) {
                            colSideColors <- rbind(colSideColorsSave,
                                                   colSideColors)
                        }
                        order2 <- order(colSideColors)
                    }
                    order1 <- order(colnames(egenefit))
                    egenefit <- egenefit[nrow(egenefit):1, ]
                    if (nrow(egenefit) == 2) {
                        names.backup <- rownames(egenefit)[1]
                        egenefit <-
                            rbind(egenefit[seq_len((nrow(egenefit)-1)), ], 0,
                                  egenefit[nrow(egenefit), ])
                        rownames(egenefit) <- seq_len(nrow(egenefit))
                        rownames(egenefit)[1] <- names.backup
                    } else {
                        egenefit <-
                            rbind(egenefit[seq_len((nrow(egenefit)-1)), ], 0,
                                  egenefit[nrow(egenefit), ])
                    }
                    rownames(egenefit)[nrow(egenefit)] <-
                        rownames(check.model)[Sgene]
                    rownames(egenefit)[(nrow(egenefit)-1)] <-
                        "---"
                    if (parameters$scoring[1] == 0) {
                        col.sums <- colMedians(abs(egenefit))
                        sig.mismatch <- which(col.sums >= 2)
                        sig.mismatch <-
                            sig.mismatch[-which(egenefit[nrow(egenefit), ]
                                                != 0)]
                        get.cols <-
                            unique(c(sig.mismatch,
                                     which(egenefit[nrow(egenefit), ]
                                           != 0)))
                        egenefit <- egenefit[, get.cols]
                    } else {
                        get.cols <- seq_len(ncol(egenefit))
                    }
                    mainlab <-
                        paste("Regulated by ", rownames(check.model)[Sgene],
                              "\n", sep = "")
                    if (csc) {
                        colSideColors <- colSideColors[get.cols]
                        if ("yellow" %in% colSideColors) {
                            mainlab <-
                                paste(mainlab,
                                      " Stimulationeffects (yellow)", sep = "")
                        }
                        if ("blue" %in% colSideColors) {
                            mainlab <-
                                paste(mainlab,
                                      " Knockdowneffects (blue)\n", sep = "")
                        }
                        if ("orange" %in% colSideColors) {
                            mainlab <-
                                paste(mainlab,
                                      " Silencingeffects Type I (orange)",
                                      sep = "")
                        }
                        if ("brown" %in% colSideColors) {
                            mainlab <-
                                paste(mainlab,
                                      " Silencingeffects Type II (brown)",
                                      sep = "")
                        }
                    } else {
                        colSideColors <- NULL
                    }
                    if (order %in% "rank") {
                        Rowv <- FALSE
                        egenefit_genes <-
                            egenefit[seq_len((nrow(egenefit)-2)), ]
                        namereset <- FALSE
                        if (!is.matrix(egenefit_genes)) {
                            egenefit_genes <- t(as.matrix(egenefit_genes))
                            namereset <- TRUE
                        } else {
                            geneorder <-
                                rownames(EtoS)[
                                    which(EtoS[, 2] == Sgene)[seq_len(Egenes)]]
                            egenefit_genes <-
                                egenefit_genes[
                                    order(match(rownames(egenefit_genes),
                                                geneorder),
                                          decreasing = TRUE), ]
                        }
                        egenefit <-
                            rbind(egenefit_genes,
                                  egenefit[(nrow(egenefit)-1):nrow(egenefit), ])
                        if (namereset) {
                            rownames(egenefit)[1] <-
                                names(which(EtoS[, 2] == Sgene))
                        }
                    }
                    if (order %in% "names") {
                        Rowv <- FALSE
                        tmp <- egenefit[seq_len((nrow(egenefit)-2)), ]
                        tmp <- tmp[order(rownames(tmp)), ]
                        rownames(egenefit)[seq_len((nrow(egenefit)-2))] <-
                            rownames(tmp)
                        egenefit[seq_len((nrow(egenefit)-2)), ] <- tmp
                    }
                    if (!is.null(dim(egenefit)) & plot) {
                        if (Rowv & nrow(egenefit) > 3) {
                            Rowv <- FALSE
                            d <- dist(egenefit[-c(nrow(egenefit)-1,
                                                  nrow(egenefit)), ])
                            hc <- hclust(d)
                            egenefit <- egenefit[c(hc$order, nrow(egenefit)-1,
                                                   nrow(egenefit)), ]
                        } else {
                            Rowv = FALSE
                        }
                        if (Colv) {
                            clusterdata <- egenefit
                            clusterdata[nrow(egenefit), ] <- 0
                            clusterdata <- NULL
                        }
                        clusterdata <- NULL
                        if ("bio" %in% colnames) {
                            colnames(egenefit) <-
                                gene2protein(
                                    myCN2bioCN(colnames(egenefit),
                                               colnames(CNOlist@stimuli),
                                               colnames(CNOlist@inhibitors)))
                        }
                        print(HeatmapOP(egenefit, main = mainlab, xrot = xrot,
                                        breaks = real.breaks, coln = 11,
                                        Colv = Colv, Rowv = Rowv,
                                        colSideColors = colSideColors,
                                        dendrogram = dendrogram, col = col,
                                        clusterx = clusterdata, ...))
                    } else {
                        print("one effect is not a matrix")
                    }
                } else {
                    print("min equals max in data matrix")
                }
            } else {
                if (csc) {
                    colSideColors <- rep("black", ncol(egenefit))
                    for (i in seq_len(length(colnames(egenefit)))) {
                        names <- unlist(strsplit(colnames(egenefit)[i], "_"))
                        if (length(names) > 1) {
                            if (names[1] %in% colnames(CNOlist@inhibitors) &
                                names[length(names)] %in%
                                colnames(CNOlist@stimuli)) {
                                colSideColors[i] <- "brown"
                            }
                            if (names[1] %in% colnames(CNOlist@stimuli) &
                                names[length(names)] %in%
                                colnames(CNOlist@inhibitors)) {
                                colSideColors[i] <- "orange"
                            }
                            if (names[1] %in% "Ctrl" & names[length(names)] %in%
                                colnames(CNOlist@stimuli)) {
                                colSideColors[i] <- "yellow"
                            }
                            if (names[1] %in% "Ctrl" & names[length(names)] %in%
                                colnames(CNOlist@inhibitors)) {
                                colSideColors[i] <- "blue"
                            }
                        } else {
                            if (names %in% colnames(CNOlist@inhibitors)) {
                                colSideColors[i] <- "blue"
                            } else {
                                colSideColors[i] <- "yellow"
                            }
                        }
                    }
                    if (!is.null(colSideColorsSave)) {
                        colSideColors <- rbind(colSideColorsSave,
                                               colSideColors)
                    }
                    order2 <- order(colSideColors)
                }
                order1 <- order(colnames(egenefit))
                egenefit <- egenefit[nrow(egenefit):1, ]
                if (nrow(egenefit) == 2) {
                    names.backup <- rownames(egenefit)[1]
                    egenefit <-
                        rbind(egenefit[seq_len((nrow(egenefit)-1)), ], 0,
                              egenefit[nrow(egenefit), ])
                    rownames(egenefit) <- seq_len(nrow(egenefit))
                    rownames(egenefit)[1] <- names.backup
                } else {
                    egenefit <-
                        rbind(egenefit[seq_len((nrow(egenefit)-1)), ], 0,
                              egenefit[nrow(egenefit), ])
                }
                rownames(egenefit)[nrow(egenefit)] <-
                    rownames(check.model)[Sgene]
                rownames(egenefit)[(nrow(egenefit)-1)] <- "---"
                if (parameters$scoring[1] == 0) {
                    col.sums <- colMedians(abs(egenefit))
                    sig.mismatch <- which(col.sums >= 2)
                    sig.mismatch <-
                        sig.mismatch[-which(egenefit[nrow(egenefit), ] != 0)]
                    get.cols <-
                        unique(c(sig.mismatch,
                                 which(egenefit[nrow(egenefit), ] != 0)))
                    egenefit <- egenefit[, get.cols]
                } else {
                    get.cols <- seq_len(ncol(egenefit))
                }
                mainlab <-
                    paste("Regulated by ", rownames(check.model)[Sgene],
                          "\n", sep = "")
                if (csc) {
                    colSideColors <- colSideColors[get.cols]
                    if ("yellow" %in% colSideColors) {
                        mainlab <- paste(mainlab,
                                         " Stimulationeffects (yellow)",
                                         sep = "")
                    }
                    if ("blue" %in% colSideColors) {
                        mainlab <- paste(mainlab,
                                         " Knockdowneffects (blue)\n",
                                         sep = "")
                    }
                    if ("orange" %in% colSideColors) {
                        mainlab <- paste(mainlab,
                                         " Silencingeffects Type I (orange)",
                                         sep = "")
                    }
                    if ("brown" %in% colSideColors) {
                        mainlab <- paste(mainlab,
                                         " Silencingeffects Type II (brown)",
                                         sep = "")
                    }
                } else {
                    colSideColors <- NULL
                }
                mainlab <- paste("Regulated by ", rownames(check.model)[Sgene],
                                 "\n", sep = "")
                egenefit[nrow(egenefit), ] <-
                    egenefit[nrow(egenefit), ]*max(abs(egenefit))
                if (order %in% "rank") {
                    Rowv <- FALSE
                    egenefit_genes <- egenefit[seq_len((nrow(egenefit)-2)), ]
                    namereset <- FALSE
                    if (!is.matrix(egenefit_genes)) {
                        egenefit_genes <- t(as.matrix(egenefit_genes))
                        namereset <- TRUE
                    } else {
                        geneorder <-
                            rownames(EtoS)[which(EtoS[, 2] ==
                                                 Sgene)[seq_len(Egenes)]]
                        egenefit_genes <-
                            egenefit_genes[order(match(rownames(egenefit_genes),
                                                       geneorder),
                                                 decreasing = TRUE), ]
                    }
                    egenefit <-
                        rbind(egenefit_genes,
                              egenefit[(nrow(egenefit)-1):nrow(egenefit), ])
                    if (namereset) {
                        rownames(egenefit)[1] <-
                            names(which(EtoS[, 2] == Sgene))
                    }
                }
                if (order %in% "names") {
                    Rowv <- FALSE
                    tmp <- egenefit[seq_len((nrow(egenefit)-2)), ]
                    tmp <- tmp[sort(rownames(tmp)), ]
                    rownames(egenefit)[seq_len((nrow(egenefit)-2))] <-
                        rownames(tmp)
                    egenefit[seq_len((nrow(egenefit)-2)), ] <- tmp
                }
                if (Rowv & nrow(egenefit) > 3) {
                    Rowv <- FALSE
                    d <- dist(egenefit[-c(nrow(egenefit)-1, nrow(egenefit)), ])
                    hc <- hclust(d)
                    egenefit <-
                        egenefit[c(hc$order, nrow(egenefit)-1,
                                   nrow(egenefit)), ]
                } else {
                    Rowv = FALSE
                }
                if (Colv) {
                    clusterdata <- egenefit
                    clusterdata[nrow(egenefit), ] <- 0
                    clusterdata <- NULL
                }
                clusterdata <- NULL
                low <-
                    sum(egenefit[nrow(egenefit), ] ==
                        min(egenefit[nrow(egenefit), ]))
                high <-
                    sum(egenefit[nrow(egenefit), ] ==
                        max(egenefit[nrow(egenefit), ]))
                egenefit2 <- egenefit
                egenefit2[which(rownames(egenefit2) %in%
                                rownames(EtoS)[which(EtoS[, 3] == -1)]), ] <-
                    egenefit2[which(rownames(egenefit2) %in%
                                    rownames(EtoS)[which(EtoS[, 3] ==
                                                         -1)]), ]*(-1)
                egenefit2 <- t(apply(egenefit2, 1, rank))
                high2 <- which(egenefit2 > ncol(egenefit2)-high)
                low2 <- which(egenefit2 < low)
                mid <- which(egenefit2 >= low & egenefit2 <=
                             ncol(egenefit2)-high)
                egenefit2[high2] <- 1
                egenefit2[low2] <- -1
                egenefit2[mid] <- 0
                if (csc) {
                    colSideColors <-
                        colSideColors[order(egenefit2[nrow(egenefit2), ])]
                }
                egenefit <- egenefit[, order(egenefit2[nrow(egenefit2), ])]
                egenefit2 <- egenefit2[, order(egenefit2[nrow(egenefit2), ])]
                show.ranks <- FALSE
                if (show.ranks) {
                    egenefit <- t(apply(egenefit, 1, rank))
                    breaks <-
                        c(0, sort(unique(egenefit[nrow(egenefit), ])),
                          ncol(egenefit)+1)
                    print(breaks)
                }
                if (ranks) {
                    if (is.null(breaks)) {
                        breaks <- c(-2,-0.5,0.5,2)
                    }
                    if ("bio" %in% colnames) {
                        colnames(egenefit2) <-
                            gene2protein(
                                myCN2bioCN(colnames(egenefit2),
                                           colnames(CNOlist@stimuli),
                                           colnames(CNOlist@inhibitors)))
                    }
                    print(HeatmapOP(egenefit2, main = mainlab, xrot = xrot,
                                    breaks = breaks, coln = 11, Colv = Colv,
                                    Rowv = Rowv, colSideColors = colSideColors,
                                    dendrogram = dendrogram,
                                    colSideColorsPos = "top", col = col,
                                    clusterx = clusterdata, ...))
                } else {
                    if (is.null(breaks)) {
                        breaks <- seq(-1,1,0.1)
                    }
                    if ("bio" %in% colnames) {
                        colnames(egenefit) <-
                            gene2protein(
                                myCN2bioCN(colnames(egenefit),
                                           colnames(CNOlist@stimuli),
                                           colnames(CNOlist@inhibitors)))
                    }
                    print(HeatmapOP(egenefit, main = mainlab, xrot = xrot,
                                    breaks = breaks, coln = 11, Colv = Colv,
                                    Rowv = Rowv, colSideColors = colSideColors,
                                    dendrogram = dendrogram,
                                    colSideColorsPos = "top", col = col,
                                    clusterx = clusterdata, ...))
                }
            }
        } else {
            mainlab <- paste("Regulated by ", rownames(check.model)[Sgene],
                             "\n", sep = "")
            if (plot) {
                print(HeatmapOP(errorMat(), main = mainlab, sub = "",
                                Colv = FALSE, Rowv = FALSE, col = "RdBu",
                                coln = 12, breaks = seq(0,8,0.1)), ...)
            }
            print("min equals max in data matrix")
            bad.data <- TRUE
        }
        if (sum(EtoS[, 2] == Sgene) == 0) {
            if (!bad.data) {
                mainlab <- paste("Regulated by ",
                                 rownames(check.model)[Sgene],
                                 "\n", sep = "")
                if (plot) {
                    print(HeatmapOP(errorMat(), main = mainlab, sub = "",
                                    Colv = FALSE, Rowv = FALSE, col = "RdBu",
                                    coln = 12, breaks = seq(0,8,0.1)), ...)
                }
                return(NULL)
            }
        } else {
            if (!is.na(rownames(as.matrix(EtoS[which(EtoS[, 2] ==
                                                     Sgene), ]))[1])) {
                if (sum(EtoS[, 2] == Sgene) > 1) {
                    genesInfo <- EtoS[which(EtoS[, 2] == Sgene), ]
                } else {
                    names.backup <- rownames(EtoS)[which(EtoS[, 2] == Sgene)]
                    genesInfo <- t(as.matrix(EtoS[which(EtoS[, 2] == Sgene), ]))
                    rownames(genesInfo) <- names.backup
                }
                if (affyIds == FALSE) {
                    temp <-
                        as.vector(unlist(mget(rownames(genesInfo),
                                              get(paste(affychip,
                                                        "SYMBOL", sep = "")))))
                    if (sum(is.na(temp) == TRUE) > 0) {
                        temp[is.na(temp)] <- rownames(genesInfo)[is.na(temp)]
                    }
                    rownames(genesInfo) <- temp
                }
                return(list(genesInfo = genesInfo, data = egenefit))
            }
        }
    }
#' sample normal form
#'
#' creates a random normal form or hyper-graph
#' @param vertices number of vertices
#' @param negation negation is allowed, if TRUE
#' @param max.edge.size maximal number of inputs per edge
#' @param max.edges maximal number of hyper-edges
#' @param dag acyclic graph, if TRUE
#' @author Martin Pirkl
#' @return random hyper-graph in normal form
#' @export
#' @examples
#' g <- randomDnf(10)
randomDnf <- function(vertices = 10, negation = TRUE, max.edge.size = NULL,
                      max.edges = NULL, dag = FALSE) {
    dnf <- NULL
    if (is.numeric(vertices)) {
        if (vertices < 27) {
            vertices <- LETTERS[seq_len(vertices)]
        } else {
            vertices <- paste("S", seq_len(vertices), "g", sep = "")
        }
    }
    if (is.null(max.edge.size)) {
        max.edge.size <- length(vertices) - 1
    }
    if (is.null(max.edges)) {
        max.edges <- length(vertices) - 1
    }
    for (i in seq_len(max.edges)) {
        edge.size <- sample(seq_len(max.edge.size), 1)
        output <- sample(vertices, 1)
        inputs <- NULL
        for (j in seq_len(edge.size)) {
            inputs <-
                c(inputs,
                  sample(vertices[-grep(paste(c(output, inputs),
                                              collapse = "|"), vertices)], 1))
        }
        if (negation) {
            pre <- sample(c("", "!"), edge.size, replace = TRUE)
        } else {
            pre <- rep("", edge.size)
        }
        inputs <- sort(inputs)
        dnf <-
            c(dnf,
              paste(c(paste(paste(pre, inputs, sep = ""), collapse = "+"), "=",
                      output), collapse = ""))
    }
    if (dag) {
        dnf <- removeCycles(dnf = dnf)
    }
    return(unique(dnf))
}
#' score a boolean network
#'
#' computes the score of a boolean network given the model and data
#' @param bString binary string denoting the boolean network
#' @param CNOlist CNOlist object
#' @param fc data (see ?bnem)
#' @param model network model
#' @param method see ?bnem
#' @author Martin Pirkl
#' @return numeric value (score)
#' @export
#' @examples
#' sim <- simBoolGtn()
#' scoreDnf(sim$bString, sim$CNOlist, sim$fc, sim$model)
scoreDnf <- function(bString, CNOlist, fc, model, method = "llr") {
    res <- bnem(search = "greedy", fc = fc, CNOlist = CNOlist,
                model = model, method = method, maxSteps = 0,
                initBstring = bString, draw = FALSE, verbose = FALSE)
    return(min(res$scores[[1]]))
}
#' plot simulation object
#'
#' plots the boolen network from a simulation
#' as disjunctive normal form
#' @param x bnemsim object
#' @param ... further arguments; see function mnem::plotDnf
#' @author Martin Pirkl
#' @return plot of boolean network
#' @method plot bnemsim
#' @export
#' @importFrom mnem plotDnf
#' @examples
#' sim <- simBoolGtn()
#' plot(sim)
plot.bnemsim <- function(x, ...) {
    plotDnf(x$model$reacID[as.logical(x$bString)])
}
#' plot bnem opbject
#'
#' plots the boolen network as disjunctive normal form
#' @param x bnemsim object
#' @param ... further arguments; see function mnem::plotDnf
#' @author Martin Pirkl
#' @return plot of boolean network
#' @method plot bnem
#' @export
#' @importFrom mnem plotDnf
#' @examples
#' sim <- simBoolGtn()
#' plot(sim)
plot.bnem <- function(x, ...) {
    plotDnf(x$model$reacID[as.logical(x$bString)])
}
