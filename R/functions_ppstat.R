# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2024-01-09
#
# Script Description: functions to handle ppstat data




# Prepare model -----------------------------------------------------------


#' Create interactions
#' 
#' Create a list of interactions as needed for a ppstat input filled
#' with zeroes.
#'
#' @param spp_names Species names vector 
#' @param times Times vector
#'
#' @return A list of lists.
#' Each list represents one species, and within each list there are the interactions
#' with this species so that flist[[i]][[j]] represents the effect of j -> i.
#' 
#' @export
#'
#' @examples
#' create_interactions_ppstat(spp_names = letters[1:5], times = seq(0, 10, by = 1))
create_interactions_ppstat <- function(spp_names, times) {
  
  # Create first level of list
  flist <- vector(mode = "list", length = length(spp_names))
  names(flist) <- spp_names
  
  # Create second (nested) level of list
  flist <- lapply(flist, function(x) flist)
  
  # Function to fill the elements of a lsit with zeroes
  fill_zeroes <- function(li, nzeroes) {
    lapply(li, function(x) rep(0, nzeroes))
  }
  
  res <- lapply(flist, function(x) fill_zeroes(x, length(times))) 
  return(res)
}

#' Write formula
#' 
#' Write a formula expected by ppstat for the inference model
#'
#' @param spp Vector of unique species names
#' @param startknot_spp First knot for species interactions. It is expressed
#' in the same unit as the stamps on the data to fit the model to (e.g. days).
#' @param endknot_spp Last knot for species interactions. It is expressed
#' in the same unit as the stamps on the data to fit the model to (e.g. days).
#' @param by_spp Step for the species interactions knots. It is expressed
#' in the same unit as the stamps on the data to fit the model to (e.g. days).
#' @param ord Order of the splines of species interactions.
#' @param hourcov Should the final formula include a hour covariate (TRUE or FALSE)?
#' @param ord_hour Order of the splines of the hour covariate for background rates.
#' @param startknot_hour First knot for the hour covariate spline for the background rate
#' (expressed in radians).
#' @param endknot_hour Last knot for the hour covariate spline for the background rate
#' (expressed in radians).
#' @param by_hour Step of the knots for the hour covariate spline for the background rate
#' (expressed in radians).
#' @param trunc Whether to truncate the species interactions. If it is NULL, species interactions splines will
#' not be truncated. Else, it must be a numerical vectors and the spline basis will be truncated to 
#' the interval from trunc[1] to trunc[2] (see the documentation of ppstat::bSpline).
#'
#' @return A string that can be interpreted as a formula.
#' 
#' @export
write_formula <- function(spp,
                          startknot_spp = -0.5, endknot_spp = 2.5, by_spp = 0.25,
                          ord = 3, 
                          hourcov = FALSE,
                          ord_hour = 3,
                          startknot_hour = -pi/4, endknot_hour = 2*pi + pi/4, by_hour = (3/24)*2*pi,
                          trunc = NULL) {
  
  # Species interaction splines ---
  if(!is.null(trunc)){
    if(length(trunc) != 2 | !is.numeric(trunc)) {
      stop("If it is not NULL trunc must be a numeric vector of length 2.")
    }
    trunc_formula <- paste0("trunc = c(", 
                            trunc[1],", ", trunc[2], ")")
  } else{
    trunc_formula <- "trunc = NULL"
  }
  
  knots_formula <- paste0("knots = seq(", 
                          startknot_spp, ", ",
                          endknot_spp, 
                          ", by = ", by_spp,")")
  ord_arg <- paste0("ord = ", ord)
  
  # Final spp spline
  spp_spline <- paste0("bSpline(x = ",
                       spp, ", ",
                       knots_formula, ", ",
                       ord_arg, ", ",
                       trunc_formula, ")")
  
  if (hourcov) {
    # Hour covariate background rate splines ---
    knots_formula_hour <- paste0("knots = seq(", 
                                 startknot_hour, ", ", 
                                 endknot_hour, 
                                 ", by = ", by_hour, ")")
    ord_arg_hour <- paste0("ord = ", ord_hour)
    hour_spline <- paste0("bSpline(x = hour, ",
                          knots_formula_hour, ", ",
                          ord_arg_hour, 
                          ", sym = TRUE)")
  }
  
  # Species response variable (write R vector with "c" in letters)
  spp_resp <- paste(spp, collapse = ", ")
  spp_resp <- paste0("c(", spp_resp, ")")
  
  # Final formula ---
  
  
  if (hourcov) {
    # Formula with hour covariates
    res <- paste0(spp_resp, " ~ ", hour_spline, " + ", paste(spp_spline, collapse = " + "))
  } else {
    # Formula without hour covariate
    res <- paste0(spp_resp, " ~ ", paste(spp_spline, collapse = " + "))
  }
  
  return(res)
}


# Data handling -----------------------------------------------------------


#' Get ppstat coefficients
#' 
#' Returns the coefficients an inferred ppstat model.
#'
#' @param model model object used for ppstat.
#' @param term dimension to pick (ie species index)
#' @param alpha wanted confidence level
#'
#' @return A dataframe with columns `coef_id`, `species`, `coef`, `lower`
#' and `upper`.
#' 
#' @export
get_ppstat_coeffs <- function(model, term, alpha = 0.05){
  
  # Get model formula
  mod.formula <- model@models[[term]]@formula
  
  # Get the index and then name of the response variable
  response.index <- attr(stats::terms(mod.formula), "response")
  response.spp <- all.vars(mod.formula)[response.index]
  
  # Get model coefficients
  coef <- model@models[[term]]@coefficients
  
  # Get coefficient table
  summ <- ppstat::summary(model)[[term]]$coefficients
  rnames <- names(rownames(summ))
  
  rownames(summ) <- NULL
  
  summ <- as.data.frame(summ)
  summ$coeff_id <- rnames
  
  e <- stats::qnorm(1-alpha/2) # epsilon
  
  # Get estimates
  coef <- summ$Estimate
  
  # Get upper and lower bounds
  error <- summ$`Std. Error`
  lower <- coef - e*error
  upper <- coef + e*error
  
  # Coef name
  coef_id <- summ$coeff_id
  
  # Response species
  species <- rep(response.spp, n=length(lower))
  
  res <- data.frame(coef_id, species, coef, lower, upper)
  
  return(res)
}

#' Get ppstat interactions
#'
#' Get toe interaction functions for ppstat
#' 
#' @param model ppstat model (assumed to be multivariate)
#' @param alpha confidence level for the confidence intervals (see the documentation of `ppstat::termPlot`)
#' @param trans transformation to apply (eg `exp` for a log-intensity) 
#' (see the documentation of `ppstat::termPlot`)
#'
#' @return A dataframe with columns:
#' 
#' + `x`: the time
#' + `variable`: response variable (species) 
#' + `value`: interaction value 
#' + `cf.lower`/`cf.upper`: confidence interval bounds
#' + `response`: response variable
#' 
#' @export
get_ppstat_interactions <- function(model, alpha = 0.05, trans = NULL) {
  
  res <- lapply(model@models, function(model) {
    pd <- ppstat:::getTermPlotData(model = model,
                                   alpha = alpha, trans = trans)
    pd$response <- ppstat:::response(model)
    return(pd)
  })
  res <- do.call(rbind, res)
}