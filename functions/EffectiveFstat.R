EffectiveFstat <- function(vX, mZ) {

  mZnorm <- mZ - mean(mZ)                                  # demean instruments
  mZZ    <- t(mZnorm) %*% mZnorm
  mZnorm <- mZnorm %*% t(chol(solve(mZZ)))                        # orthogonalize instruments 
  
  vCoefs <- solve(t(mZnorm) %*% mZnorm) %*% t(mZnorm) %*% vX      # ols reduced form 
  vRes   <- vX - mZnorm %*% vCoefs                                # residuals   
  
  mW     <- NeweyWestRes(vRes, mZnorm ,4)                         # variance covariance  
  
  dF     <- t(vX) %*% mZnorm %*% t(mZnorm) %*% vX / sum(diag(mW)) # F-stat    
  
  return(dF)
}
