## transRed <- function(g) { # general transitive reduction:
##   g2 <- g
##   for (i in g) {
##     wtw <- unlist(strsplit(i, "="))
##     for (j in g[grep(paste("=", wtw[2], sep = ""), g)]) {
##       wtw2 <- unlist(strsplit(j, "="))
##       if (length(grep(paste(wtw[1])))) {}
##   return(g)
## }
