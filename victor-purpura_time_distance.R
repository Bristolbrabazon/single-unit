vp.time.distance <- function(spiketrain1, spiketrain2, cost) {
  # The function calculates the "spike time" distance (Victor & Purpura 1996) for a single cost
  # 
  # spiketrain1: vector of spike times for first spike train
  # spiketrain2: vector of spike times for second spike train
  # cost: cost per unit time to move a spike
  # 
  # Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
  # Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
  # Translated to R (4.2.3) by Daniil Kiselev.
  
  n_sp_i <- length(spiketrain1)
  n_sp_j <- length(spiketrain2)
  
  if(cost == 0) {
    return(abs(n_sp_i-n_sp_j))
  } else if(cost == Inf) {
    return(n_sp_i+n_sp_j)
  }
  
  scr <- matrix(data = 0,
                nrow = n_sp_i+1,
                ncol = n_sp_j+1) # Initializes margins with cost of adding a spike
  scr[, 1] <- 0:n_sp_i
  scr[1, ] <- 0:n_sp_j
  
  if(n_sp_i && n_sp_j) {
    for(i in 2:(n_sp_i+1)) {
      for(j in 2:(n_sp_j+1)) {
        scr[i, j] <- min(scr[i-1,j]+1, scr[i,j-1]+1, scr[i-1,j-1]+cost*abs(spiketrain1[i-1]-spiketrain2[j-1]))
      }
    }
  }
  
  scr[n_sp_i+1, n_sp_j+1]
}
