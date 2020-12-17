FindUpperLower <- function(vTest, vPar) {

  vBounds <- c(0, 0)
  
  for (k in 0:length(vPar)) {
    if (vTest[k+1] != 0) {
        vBounds[1] <- vPar[k+1]
      break
    }
    if (vTest[k+1] != 0 && vTest[k+2] == 0) {
      vBounds[1] <- vPar[k+2]
      break
    }
  }

  for (k in length(vPar):2) {
    if (vTest[k] != 0) {
      vBounds[2] <- vPar[k]
      break
    }
    if (vTest[k] != 0 && vTest[k-1] == 0) {
      vBounds[2] <- vPar[k-1]
      break
    }
  }
  
  vBounds
}
