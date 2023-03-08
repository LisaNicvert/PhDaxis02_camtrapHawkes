# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2023-02-28
#
# Script Description: functions to compute sensitivity and specificity and
# to format results dataframes containing inference methods performance

compute_sensi <- function(TP, FN) {
  # Computes sensitivity
  # ### Inputs
  # TP: true positives vector
  # FN: false negatives vector
  # ### Outputs
  # The sensitivity (true positive rate)
  
  sensi <- TP/(TP+FN)
  return(sensi)
}

compute_speci <- function(TN, FP) {
  # Computes specificity
  # ### Inputs
  # TN: true negatives vector
  # FP: false positives vector
  # ### Output
  # The specificity (true negative rate)
  speci <- TN/(TN+FP)
  return(speci)
}

format_data_perf <- function(d, add_quantiles = TRUE, level = 0.05){
  # Format a dataframe measuring sensitivity and specificity
  # ### Inputs
  # dataframe d:
  #   must have columns TP, FP, TN, FN.  additional columns possible, 
  #   but additional rows will affect quantiles computation if add_quantiles is TRUE.
  # add_quantiles: should quantiles for sensitivity and specificity be computed?
  #   Quantiles are computed on values grouped by all columns (except "TP", "FP", "TN", "FN").
  # level = level of the quantiles.
  # ### Output
  # Returns a cleaned form of d:
  #   additional column sensi computed as TP/(TP+FN)
  #   additional column speci computed as TN/(TN+FP)
  #   the columns are pivoted to longer format so that
  #   sensi and speci are grouped into "value" column and
  #   "type" describes sensi or speci.
  #   if add_quantiles, quantile values of value (sensi/speci) will
  #   be added (quantile computed grouping according to all other columns 
  #   excluding "TP", "FP", "TN", "FN").
  
  dres <- d
  
  # Compute sensitivity and specificity
  dres <- dres %>% mutate(sensi = compute_sensi(TP, FN),
                          speci = compute_speci(TN, FP))
  # Pivot data
  dres <- dres %>% 
    pivot_longer(cols = c("sensi", "speci"), names_to = "type")
  
  if(add_quantiles){
    other_cols <- colnames(dres)[which(!(colnames(dres) %in% c("TP", "FP", "TN", "FN", "value")))]
    
    dres <- dres %>% 
      group_by(across(all_of(other_cols))) %>%
      mutate(median = median(value)) %>%
      mutate(qinf = quantile(value, level/2)[[1]]) %>%
      mutate(qsup = quantile(value, 1-(level/2))[[1]])
  }
  
  return(dres)
}
