checkMethod <-
function(method) {
  methods2 <- method
  methods <- c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "spearman", "pearson", "kendall", "mLL", "cp", "none")
  method <- methods[grep(paste(paste("^", method, sep = ""), collapse = "|"), methods)[1]]
  method  <- unique(c(method, methods2))
  if (!(any(c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "spearman", "pearson", "kendall", "mLL", "cp", "none") %in% method))) {
    stop("I can't let you do that, Dave. You have to pick a valid method: cosine, euclidean, maximum, manhattan, canberra, binary, minkowski, spearman, pearson, kendall, mLL, cp, none")
  }
  return(method)
}
