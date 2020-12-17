NeweyWestRes <- function(e, X, L) {
  
  N <- dim(X)[1]
  k <- dim(X)[2]
  
  if (L < 0) {
    L <- floor(4*((N/100)^(2/9)))
  }
  
  Q <- 0
  for (l in 0:L) {
    w_l  <- 1 - l / (L + 1)
    for (t in (l+1):N) {
      if (l == 0) {   
        # This calculates the S_0 portion
        Q = Q  + e[t]^2 * t(X[t, ]) %*% X[t, ]
      } else {        
        # This calculates the off-diagonal terms
        Q = Q + w_l * e[t] * e[t-1] * (t(X[t, ]) %*% X[t-l, ] + t(X[t-l, ]) %*% X[t, ])
      }
    }
  }
  Q
}
