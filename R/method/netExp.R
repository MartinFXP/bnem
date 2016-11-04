netExp <- function(n) {
  exp <- 0
  for (i in 1:n) {
    add <- n*choose(n,i)*2^i
    sub <- n*choose(n*2-i,i-1) # i*n
    exp <- exp + add - sub
  }
  return(exp)
}

netExpFull <- function(n) {
  exp <- 0
  for (i in 1:(n*2)) {
    add <- choose(n*2,i)
    exp <- exp + add
  }
  return(exp)
}

## for (i in 1:10) {
##   print(netExp(i))
##   print(netExpFull(i))
##   print(netExpFull(i) - netExp(i))
##   print(netExp(i)/netExpFull(i))
## }
