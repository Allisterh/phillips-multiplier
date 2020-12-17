PhillipsMultiplier <- function(vParB, iH, vPi, vU, mZ, mX, CI) {
  
  vBeta       <- rep(0, iH)  # storage 2SLS point PM 
  vBetauc     <- rep(0, iH)  # storage OLS point PM 
  vBetaPi     <- rep(0, iH)  # storage irf pi   
  vBetaUr     <- rep(0, iH)  # storage irf ur        
  vLowerAR    <- rep(0, iH)  # lower AR-bound  
  vUpperAR    <- rep(0, iH)  # upper AR-bound 
  vLowerAR2s  <- rep(0, iH)  # lower 2SLS-bound  
  vUpperAR2s  <- rep(0, iH)  # upper 2SLS-bound 
  vLowerARuc  <- rep(0, iH)  # lower unconditional multiplier 
  vUpperARuc  <- rep(0, iH)  # upper unconditional multiplier 
  vLowerPI    <- rep(0, iH)  # lower bound AR pi irf
  vUpperPI    <- rep(0, iH)  # upper bound AR pi irf 
  vLowerUR    <- rep(0, iH)  # lower bound AR ur irf
  vUpperUR    <- rep(0, iH)  # upper bound AR ur irf
  vFstat      <- rep(0, iH)  # F-stats 
  
  mWald       <- matrix(0, length(vParB), iH)   # store Wald
  mTest       <- matrix(0, length(vParB), iH)   # store rejection 
  
  for (h in 1:iH) {
    
    vPiC <- vPi   
    vUC  <-  vU
    
    # construct cumulative pi and u
    if (h > 1) {
      for (j in 2:h) {
        vPiC <- vPiC + LagMatrix(vPi, -(j-1))
        vUC  <- vUC  + LagMatrix(vU, -(j-1))
      }
    }
    
    mObs <-  cbind(vPiC, vUC, mZ, mX) %>% na.omit()
    
    iTc   <- dim(mObs)[1]
    iCols <- dim(mObs)[2]
    
    mXC  <- cbind(rep(1, iTc), mObs[,4:iCols])
    mMx  <- diag(iTc) - mXC %*% solve(t(mXC) %*% mXC) %*% t(mXC)
    vPiC <- mMx %*% mObs[, 1]
    vUC  <- mMx %*% mObs[, 2]
    mZC  <- mObs[, 3] %>% as.matrix()  
    
    vCoefs   = solve( t(vUC) %*% mZC %*% solve(t(mZC) %*% mZC) %*% t(mZC) %*% vUC) %*% t(vUC) %*% mZC %*% solve(t(mZC) %*% mZC) %*% t(mZC) %*% vPiC 
    vBeta[h] = vCoefs[1]
    
    vRes          = vPiC -  vBeta[h] * vUC
    mV            = NeweyWestRes(vRes, mZC , -1)
    mPz           = mZC %*% solve(t(mZC) %*% mZC) %*% t(mZC)
    mAvar         = solve(t(vUC) %*% mPz %*% vUC) %*% t(vUC) %*% mZC %*% solve(t(mZC) %*% mZC) %*% mV * solve(t(mZC) %*% mZC) %*% t(mZC) %*% vUC %*% solve(t(vUC) %*% mPz %*% vUC)
    vLowerAR2s[h] = vBeta[h] - qnorm(CI) * sqrt(mAvar)
    vUpperAR2s[h] = vBeta[h] + qnorm(CI) * sqrt(mAvar)
    
    # AR type CS 
    for (k in 1:length(vParB)) {
      vRes      = vPiC - vParB[k] * vUC
      vGamma    = t(mZC) %*% vRes
      vRes1     = vRes - mZC %*%  solve(t(mZC) %*% mZC) %*% vGamma       
      mV        = NeweyWestRes(vRes1, mZC , -1)
      mWald[k,h] = t(vGamma) %*% solve(mV) %*% vGamma
      if (mWald[k,h] < qchisq(CI, 1)) {
        mTest[k,h] <- 1 
      }
    }
    
    # find upper and lower confidence bounds AR statistic  
    vBounds     <- FindUpperLower(mTest[, h], vParB)
    vLowerAR[h] <- vBounds[1]
    vUpperAR[h] <- vBounds[2]    
    
    # compute effective F-stat 
    vFstat[h] <- EffectiveFstat(vUC, mZC)
    
    # Compute IRF PI 
    vBetaPi[h] <- solve(t(mZC) %*% mZC) %*% t(mZC) %*% vPiC
    vRes       <- vPiC - mZC *  vBetaPi[h]       
    mV         <- NeweyWestRes(vRes, mZC ,4)
    mAvar      <- solve(t(mZC) %*% mZC) %*% mV %*% solve(t(mZC) %*% mZC)
    vLowerPI[h]<- vBetaPi[h] - qnorm(CI) * sqrt(mAvar)
    vUpperPI[h]<- vBetaPi[h] + qnorm(CI) * sqrt(mAvar)
    
    # Compute IRF UR
    vBetaUr[h] <- solve(t(mZC) %*% mZC) %*% t(mZC) %*% vUC
    vRes       <- vUC - mZC *  vBetaUr[h]       
    mV         <- NeweyWestRes(vRes, mZC ,4)
    mAvar      <- solve(t(mZC) %*% mZC) %*% mV %*% solve(t(mZC) %*% mZC)
    vLowerUR[h]<- vBetaUr[h] - qnorm(CI) * sqrt(mAvar)
    vUpperUR[h]<- vBetaUr[h] + qnorm(CI) * sqrt(mAvar)
    
    # Compute unconditional PM 
    vCoefs        <- solve(t(vUC) %*% vUC) %*% t(vUC) %*% vPiC
    vBetauc[h]    <- vCoefs[1]
    vRes          <- vPiC - vBetauc[h] * vUC
    mV            <- NeweyWestRes(vRes, vUC , -1)
    Avar          <- solve(t(vUC) %*% vUC) %*% mV %*% solve(t(vUC) %*% vUC)
    vLowerARuc[h] <- vBetauc[h] - qnorm(CI) * sqrt(Avar)
    vUpperARuc[h] <- vBetauc[h] + qnorm(CI) * sqrt(Avar)   
  }

  h <- 1:iH
  
  cbind(
    h,
    vBeta,
    vUpperAR,
    vLowerAR,
    vBetauc,
    vUpperARuc,
    vLowerARuc,
    vBetaPi,
    vUpperPI,
    vLowerPI,
    vBetaUr,
    vUpperUR,
    vLowerUR,
    vFstat
  )
  
}
