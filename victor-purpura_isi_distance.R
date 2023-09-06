vp.isi.distance <- function(spiketrain1, spiketrain2, cost, int_length) {
  # Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
  # Translated to Matlab by Christina Behrend from FORTRAN code by Jonathan Victor.
  # Matlab code optimized for faster processing by Jim Hokanson
  # Translated to R (4.2.3) by Daniil Kiselev
  # 
  # calculates distance between two spike trains in the spike interval metric
  # by a continuum modification of the sellers algorithm
  # 
  # end conditions: the first and last ISI are expanded as needed to minimize
  # the total cost
  # 
  # input variables
  # spiketrain1: vector spike times for first spike train
  # spiketrain2: vector spike times for second spike train
  # cost: cost per unit time to move a spike
  # int_length: the length of the entire window (in datapoints)
  
  n_sp_i <- length(spiketrain1)
  n_sp_j <- length(spiketrain2)
  
  n_i <- n_sp_i+1
  n_j <- n_sp_j+1
  scr <- matrix(data = 0,
                nrow = n_i+1,
                ncol = n_j+1)
  
  if(cost == 0) {
    return(abs(n_i-n_j))
  }
  
  scr[, 1] <- 0:n_i
  scr[1, ] <- 0:n_j
  
  sti_diff <- diff(spiketrain1)
  stj_diff <- diff(spiketrain2)
  
  for(i in 1:n_i) {
    if((i > 1)&&(i < n_i)) {
      d_i <- sti_diff[i-1]
    } else if((i == 1)&&(i == n_i)) {
      d_i <- int_length
    } else if((i == 1)&&(i < n_i)) {
      d_i <- spiketrain1[i]
    } else {
      d_i <- int_length-spiketrain1[i-1]
    }
    
    i_end <- ((i == 1)||(i == n_i))
    
    if(n_j == 1) {
      d_j <- int_length
    } else {
      d_j <- spiketrain2[1]
    }
    
    if(i_end) {
      dist <- 0
    } else {
      dist <- max(0, d_j-d_i)
    }
    
    scr[i+1, 2] <- min(scr[i, 2]+1, scr[i+1, 1]+1, scr[i, 1]+cost*dist)
    
    for(j in 2:(n_j-1)) {
      d_j <- stj_diff[j-1]
      
      if(i_end) {
        dist <- max(0, d_i-d_j)
      } else {
        dist <- abs(d_i-d_j)
      }
      
      tmp <- min(scr[i, j+1]+1, scr[i+1, j]+1)
      scr[i+1, j+1] <- min(tmp, scr[i, j]+cost*dist)
    }
    
    if(n_j == 1) {
      d_j <- int_length
    } else {
      d_j <- int_length-spiketrain2[n_j-1]
    }
    
    if(i_end) {
      dist <- 0
    } else {
      dist <- max(0, d_j-d_i)
    }
    
    scr[i+1, n_j+1] <- min(scr[i, n_j+1]+1, scr[i+1, n_j]+1, scr[i, n_j]+cost*dist)
  }
  
  scr[n_i+1, n_j+1]
}
