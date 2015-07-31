preClingo <- function(NEMlist, model, stimuli, inhibitors, path = ".\\") {
  batches <- NULL
  runs <- NULL
  
  vertices <- c(stimuli, inhibitors, "Dump1", "Dump2")
  stimuli <- stimuli
  inhibitors <- inhibitors
  readouts <- c(inhibitors,"Dump1")
  edges <- model$reacID[-grep("\\+", model$reacID)]
  exps <-makeDesign(NEMlist$exprs, stimuli, inhibitors, batches = c(batches, runs))
  exps <- exps[, c(stimuli, inhibitors)]
  obs <- NEMlist$exprs
  obs.disc <- disc(simpleNorm(obs))
  fcs <- NEMlist$fc
  fcs.disc <- disc(fcs)
  NEMlist$norm <- obs.disc

  write("% start",  file = paste(path, "input.lp", sep = ""), append = FALSE)
  for (i in vertices) {
    write(paste("vertex(", tolower(i), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  for (i in stimuli) {
    write(paste("stimuli(", tolower(i), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  for (i in inhibitors) {
    write(paste("inhibitor(", tolower(i), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  for (i in readouts) {
    write(paste("readout(", tolower(i), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  for (i in edges) {
    nodes <- strsplit(i, "=")
    nodes <- unlist(nodes)
    if (length(grep("!", nodes[1])) == 1) {
      type = -1
      nodes[1] <- gsub("!", "", nodes[1])
    } else {
      type = 1
    }
    write(paste("edge(", tolower(nodes[1]), ",", tolower(nodes[2]), ",", type, ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  for (i in 1:nrow(exps)) {
    for (j in 1:max(sum(exps[i, ] == 1),1)) {
      if (sum(exps[i, ] == 1) == 0) {
        for (k in stimuli) {
          write(paste("exp(", i, ",", tolower(k), ",", 0, ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
        }
        next()
      }
      for (k in stimuli) {
        if (exps[i, k] == 1) {
          write(paste("exp(", i, ",", tolower(k), ",", 1, ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
        } else {
          write(paste("exp(", i, ",", tolower(k), ",", 0, ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
        }
      }
      node <- tolower(colnames(exps)[which(exps[i, ] == 1)[j]])
      if (node %in% tolower(inhibitors)) {
        type <- 0
        write(paste("exp(", i, ",", node, ",", type, ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
      }
    }
  }
  write(paste("nexp(", ncol(obs.disc), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  for (i in 1:nrow(NEMlist$exprs)) {
    write(paste("egene(", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(NEMlist$exprs)[i]))))), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  }
  write(paste("negenes(", nrow(obs.disc), ").", sep = ""), file = paste(path, "input.lp", sep = ""), append = TRUE)
  
  write("% start", file = paste(path, "data.lp", sep = ""), append = FALSE)
  for (i in 1:nrow(obs.disc)) {
    for (j in 1:ncol(obs.disc)) {
      if (obs.disc[i, j] == 1) {
        write(paste("obs(", j, ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rownames(obs.disc)[i]))))), ").", sep = ""), file = paste(path, "data.lp", sep = ""), append = TRUE)
      }
    }
  }

  write("% start", file = paste(path, "dataOld.lp", sep = ""), append = FALSE)
  for (i in 1:nrow(obs.disc)) {
    write(paste("obs(", 1:ncol(obs.disc), ",", gsub("\\:", "_", gsub("\\*", "_", gsub("\\-", "_", gsub("\\+", "_", tolower(rep(rownames(obs.disc)[i], ncol(obs.disc))))))), ",", obs.disc[i, ], ").", sep = ""), file = paste(path, "dataOld.lp", sep = ""), append = TRUE)
  }

  write("% start", file = paste(path, "dataFc.lp", sep = ""), append = FALSE)
  NEMlist$fc <- computeFcAsp(CNOlist, NEMlist$exprs, file = paste(path, "dataFc.lp", sep = ""))
  print("data successfully converted; open terminal and type:")
  print(paste("clingo --opt-all 0 ", path, "input.lp ", path, "data.lp ", path, "program.lp", sep = ""))
}

#./clingo --opt-all 0 input.lp data.lp booNemAbs.lp
# clingo_64 --restart-on-model 0 input.lp dataFc.lp booNemFc.lp
# clingo_64 --restart-on-model 0 input.lp data.lp booNemAbs.lp
