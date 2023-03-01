# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-12-05
#
# Script Description: functions to compute intervals between events

compute_intervals <- function(sp_from, sp_to) {
  # Compute the inter-event times. Function adapted from
  # Murphy et al 2021 (http://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13548)
  # ### Input
  # sp_from: species that is observed first
  # sp_to: following species
  # ### Output
  # Returns the median of the time interval between the 2 
  # species occurrences
  
  # Get observed inter-event times
  maxrow <- max(nrow(sp_from), nrow(sp_to))
  
  if(maxrow == nrow(sp_from)) {
    maxobs <- sp_from
    minobs <- sp_to
  } else {
    maxobs <- sp_to
    minobs <- sp_from
  }
  
  maxobs <- as.data.frame(maxobs)
  minobs <- as.data.frame(minobs)
  
  fromto <- rep(NA, maxrow)
  
  for (i in 1:maxrow) {
    # subsetting the detections of spp with the least
    # detections seen where the more abd spp was seen
    min.temp <- minobs[minobs$camera == maxobs[i,"camera"],]
    if (nrow(min.temp)>0) {
      if(maxrow == nrow(sp_from)) { # if maxobs is the first spp
        hold <- min.temp$stamp - maxobs[i,"stamp"]
      } else { # if minobs is the first species
        hold <- maxobs[i,"stamp"] - min.temp$stamp
      }
      
      fromto[i]  <- suppressWarnings(min(hold[hold>0])) # will be NA if to is never seen after from
    } else {fromto[i] = "Nodata"}
  }
  
  # Make the characters numbers
  fromto <- suppressWarnings(as.numeric(fromto))
  
  # Make infinite values NAs
  is.na(fromto) <- sapply(fromto, is.infinite)
  
  return(fromto)
}

compute_intervals_permute <- function(sp_from, sp_to) {
  # Compute the permuted inter-event times. Function adapted from
  # Murphy et al 2021 (http://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13548)
  # ### Input
  # sp_from: species that is observed first
  # sp_to: following species
  # ### Output
  # Returns the median of the permuted time interval between the 2 
  # species occurrences
  
  # Get observed inter-event times
  maxrow <- max(nrow(sp_from), nrow(sp_to))
  
  if(maxrow == nrow(sp_from)) {
    maxobs <- sp_from
    minobs <- sp_to
  } else {
    maxobs <- sp_to
    minobs <- sp_from
  }
  
  maxobs <- as.data.frame(maxobs)
  minobs <- as.data.frame(minobs)
  
  times.sims <- rep(NA, maxrow)
  #subsetting locations based on stations that are available at the same time as the detection
  for (i in 1:maxrow){
    loc.temp <- subset(maxobs,
                       maxobs$stamp[i] >= maxobs$Begin & maxobs$stamp[i] <= maxobs$End)
    ##shuffling among the subsetted stations
    sim.max <- sample(loc.temp$camera, 
                      length(loc.temp$camera), 
                      replace = FALSE)
    
    ###choosing a random subsetted station
    min.temp <- minobs[minobs$camera == sim.max[1],]
    if (nrow(min.temp) > 0) {
      # minimum time to capture
      if(maxrow == nrow(sp_from)) { # if maxobs is the first spp
        sims.hold <- min.temp$stamp - maxobs$stamp[i]
      } else { # if minobs is the first species
        sims.hold <- maxobs$stamp[i] - min.temp$stamp
      }
      times.sims[i] <- suppressWarnings(min(sims.hold[sims.hold > 0]))
    } else {times.sims[i] = "Nodata"}
  }
  times.sims <- suppressWarnings(as.numeric(times.sims))
  times.sims[is.infinite(times.sims)] <- NA
  
  return(times.sims)
}

compute_TN_TP <- function(pval, true, alpha) {
  # Compute True positive and true negative rates from p-values.
  # ### Inputs
  # pval: a vector of p-values
  # true: the ground truth vector (must be in the same order as pval)
  # alpha: the significance threshold (defaults to 0.05)
  # ### Output
  # A named vector TP, FP, TN, FN with the values.
  
  
  # Transform p-values to 0/1 based on comparison with alpha
  inferred <- as.numeric(pval <= alpha)
  # 0 = not inferred
  # 1 = inferred
  
  # -- Positives
  Pos <- inferred[true != 0]
  
  # True Positives
  TP <- length(Pos[Pos != 0])
  
  # False negatives
  FN <- length(Pos[Pos == 0])
  
  # -- negatives
  Neg <- inferred[true == 0]
  
  # True negatives
  TN <- length(Neg[Neg == 0])
  
  # False positives
  FP <- length(Neg[Neg != 0]) 
  
  res <- c("TP" = TP, "FP" = FP, "TN" = TN, "FN" = FN)
  return(res)
}
