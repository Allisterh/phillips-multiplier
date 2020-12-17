LagMatrix <- function (x, lag) 
{
  n <- length(x)
  k <- length(lag)
  mlg <- max(c(0, lag[lag > 0]))
  mld <- max(abs(c(0, lag[lag < 0])))
  lmat <- array(NA, c(n + mlg + mld, k))
  for (i in 1:k) {
    lmat[(1 + lag[i] + mld):(n + lag[i] + mld), i] <- x
  }
  lmat <- lmat[(mld + 1):(mld + n), , drop = FALSE]
  return(lmat)
}