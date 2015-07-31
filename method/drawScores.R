drawScores <- function(CNOresult) {
  gens <- CNOresult$results[, 1]
  stallGens <- as.numeric(unlist(strsplit(gens, " / ")))
  stallGens <- stallGens[(1:(length(stallGens)/2))*2 - 1]
  getInfo <- CNOresult$results[which(stallGens == 0), ]
  bestScore <- unlist(strsplit(getInfo[, 2], " "))
  bestScore <- as.numeric(gsub("\\(", "", gsub("\\)", "", bestScore)))
  bestScore <- bestScore[(1:(length(bestScore)/2))*2 - 1]
  avgScore <- unlist(strsplit(getInfo[, 2], " "))
  avgScore <- as.numeric(gsub("\\)", "", gsub("\\(", "", avgScore[(1:(length(avgScore)/2))*2])))
  par(mfrow=c(3,1))
  plot(1:length(bestScore), bestScore, col = "red", type = "l", main = paste("Score Improvement", sep = ""), xlab = "Generation", ylab = "Score", ylim = c(min(bestScore), max(c(bestScore, avgScore))), xaxt='n')
  axis(1, at = 1:length(bestScore), labels = which(stallGens == 0))
  lines(1:length(bestScore), avgScore, col = "blue", type = "l")
  legend(length(bestScore), max(c(bestScore, avgScore)), legend = c("Best Score", "Average Score"), fill = c("red", "blue"), xjust = 1)
  if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 200) {
    overallTime <- paste(round(sum(as.numeric(CNOresult$results[, "Iter_time"]))/60, 3), " minutes", sep = "")
  }
  if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 18000) {
    overallTime <- paste(round(sum(as.numeric(CNOresult$results[, "Iter_time"]))/60/60, 4), " hours", sep = "")
  }
  if (sum(as.numeric(CNOresult$results[, "Iter_time"])) > 864000) {
    overallTime <- paste(round(sum(as.numeric(CNOresult$results[, "Iter_time"]))/60/60/24, 5), " days", sep = "")
  }
  if (sum(as.numeric(CNOresult$results[, "Iter_time"])) <= 200 ) {
    overallTime <- paste(round(sum(as.numeric(CNOresult$results[, "Iter_time"])), 2), " seconds", sep = "")
  }
  plot(as.numeric(CNOresult$results[, "Iter_time"]), xlab = "Generation", ylab = "Iteration Time in seconds", main = paste("Iteration Time over all Generations (time overall: ", overallTime, ")", sep = ""), type = "l")
  #lines(lowess(as.numeric(CNOresult$results[, "Iter_time"])), col = "red")
  bStringDist <- numeric()
  for (i in 1:(nrow(getInfo) - 1)) {
    bStringDist <- c(bStringDist, dist(rbind(unlist(strsplit(getInfo[i, 3], ", ")), unlist(strsplit(getInfo[i+1, 3], ", ")))))
  }
  bStringDist <- c(0, bStringDist)
  plot(1:length(bStringDist), bStringDist, type = "l", main = paste("Best Strings Distance (average distance: ", round(mean(bStringDist), 2), ")", sep = ""), xlab = "Generation", ylab = "Euclidean Distance", xaxt='n')
  #mtext("Edge Difference", 4, line = 2)
  axis(1, at = 1:length(bestScore), labels = which(stallGens == 0))
  diffs <- numeric()
  for (i in 1:length(unlist(strsplit(getInfo[1, 3], ", ")))) {
    diffs <- c(diffs, dist(rbind(rep(0, i), rep(1, i))))
  }
  #abline(h = diffs, lty = 3)
  #abline(v = bStringDist, lty = 3)
  axis(4, at = diffs, labels = 1:length(diffs))
  #as.numeric(unlist(strsplit(temp, ", "))) # change taht°!!!!!!!!°°°°°°°°°°°°°°°°°°°°°°°°°°°¹²³¼½¬{[]}\¸¬@ł€¶ŧ←↓→øþ¨~`^˝łĸŋđðſæ|»«¢„“”µ·…
}
