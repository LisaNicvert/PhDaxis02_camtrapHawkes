# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2023-02-28
#
# Script Description: functions to compute sensitivity and specificity and
# to format results dataframes containing inference methods performance

#' Computes sensitivity
#' 
#' Computes sensitivity from true positives and false negatives
#'
#' @param TP true positives vector
#' @param FN false negatives vector
#'
#' @return The sensitivity (true positive rate)
#' @export
compute_sensi <- function(TP, FN) {
  sensi <- TP/(TP+FN)
  return(sensi)
}

#' Computes specificity
#' 
#' Computes specificity from true negatives and false positives
#'
#' @param TN true negatives vector
#' @param FP false positives vector
#'
#' @return The specificity (true negative rate)
#' @export
compute_speci <- function(TN, FP) {
  speci <- TN/(TN+FP)
  return(speci)
}

#' Format performance dataframe
#' 
#' Format a dataframe measuring sensitivity and specificity
#'
#' @param d dataframe: must have columns TP, FP, TN, FN.  
#' additional columns possible, but they will affect 
#' quantiles computation if add_quantiles is TRUE.
#' @param add_quantiles should quantiles for sensitivity and specificity be computed?
#' Quantiles are computed on values grouped by all columns (except "TP", "FP", "TN", "FN").
#' @param level level of the quantiles.
#'
#' @return Returns a cleaned form of d:
#'   additional column sensi computed as TP/(TP+FN)
#'   additional column speci computed as TN/(TN+FP)
#'   the columns are pivoted to longer format so that
#'   sensi and speci are grouped into "value" column and
#'   "type" describes sensi or speci.
#'   if add_quantiles, quantile values of value (sensi/speci) will
#'   be added (quantile computed grouping according to all other columns 
#'   excluding "TP", "FP", "TN", "FN").
#' @export
format_data_perf <- function(d, add_quantiles = TRUE, level = 0.05){
  
  dres <- d
  
  # Compute sensitivity and specificity
  dres <- dres %>% mutate(sensi = compute_sensi(TP, FN),
                          speci = compute_speci(TN, FP))
  # Pivot data
  dres <- dres %>% 
    tidyr::pivot_longer(cols = c("sensi", "speci"), names_to = "type")
  
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
