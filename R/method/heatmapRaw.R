heatmapRaw <- function(matrix,
                       main = "heatmapRaw",
                       cluster = 0,
                       dendrogram = "none",
                       col = "redgreen",
                       scale = "none",
                       ...) {
  require("gplots")
  if (cluster == 3) {
    Rowv = TRUE
    Colv = TRUE
  }
  if (cluster == 2) {
    Rowv = TRUE
    Colv = FALSE
  }
  if (cluster == 1) {
    Rowv = FALSE
    Colv = TRUE
  }
  if (cluster == 0) {
    Rowv = FALSE
    Colv = FALSE
  }
  heatmap.2(matrix,
            Colv = Colv,
            Rowv = Rowv,
            trace = "none",
            # colsep = colsep, rowsep = FALSE, sepcol = "White", sepwidth = rep(0.01, samples/t - 1),
            main = main,
            dendrogram = dendrogram,
            col = col,
            scale = scale,
            ...)
}
