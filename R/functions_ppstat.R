get_ppstat_coeffs <- function(model, term, alpha = 0.05){
  # Returns the coefficients of a ppstat model.
  # model is the model object used for ppstat.
  # term is the dimension to pick (ie species index)
  # alpha is the wanted confidence level (default = 5 %)
  
  # Get model formula
  mod.formula <- model@models[[term]]@formula
  
  # Get the index and then name of the response variable
  response.index <- attr(terms(mod.formula), "response")
  response.spp <- all.vars(mod.formula)[response.index]
  
  # Get model coefficients
  coef <- model@models[[term]]@coefficients
  
  # Get coefficient table
  summ <- summary(model@models[[term]])$coefficients
  rnames <- names(rownames(summ))
  
  rownames(summ) <- NULL
  
  summ <- as.data.frame(summ)
  summ$coeff_id <- rnames
  
  e <- qnorm(1-alpha/2) # epsilon
  
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


my_getTermPlotData <- function(model, alpha = 0.05, trans = NULL, ...) {
  # Helper function for getPlotData
  # Rewritten from getTermPlotData in ppstat package
  
  if (alpha <= 0 || alpha > 1) # Check alpha = confidence level
    stop("The 'alpha' level must be in (0,1]")
  if (isTRUE(all.equal(alpha, 1))) { # check if confidence level is one (ie no confidence interval)
    se <- FALSE
  }else{
    se <- TRUE
    q <- qnorm(1-alpha/2)
  }
  
  linearFilter <- getLinearFilter(model, se = se, nr = 400)
  
  if (se){ # if a confidence interval is needed
    moltenFilter <- reshape2::melt(linearFilter$linearFilter, id.vars = "x")
    plotData <- cbind(moltenFilter,
                      data.frame(cf.lower = moltenFilter$value - q*unlist(linearFilter$se),
                                 cf.upper = moltenFilter$value + q*unlist(linearFilter$se)))
    if (!is.null(trans))
      plotData[, c("value", "cf.lower", "cf.upper")] <- do.call(trans, 
                                                                list(plotData[, c("value", 
                                                                                  "cf.lower", 
                                                                                  "cf.upper")]))
  }else{ # no confidence interval
    plotData <- reshape2::melt(linearFilter, id.vars = "x")
    if(!is.null(trans))
      plotData$value <- do.call(trans, plotData$value)
  }
  return(plotData)
}



getPlotData <-function (model, ...){
  # Function to return data used for termPlots in ppstat
  # Draws heavily on termPlot function from ppstat package.
  # returns a dataframe with columns 
  #   x = time (delay)
  #   variable = control variables (species from)
  #   value = value of excitation function
  #   response = response variable (species to)
  .local <- function(model, alpha = 0.05, layer = geom_line(), 
                     trans = NULL, ...){
    noFilterModels <- sapply(model@models, function(m) length(m@filterTerms) == 0)
    if (all(noFilterModels)) {
      print("No filter function terms to plot.")
      return(invisible())
    }
    plotData <- lapply(model@models[!noFilterModels], function(model) {
      pd <- my_getTermPlotData(model = model, alpha = alpha, 
                               trans = trans, ...)
      pd$response <- do.call(function(...) paste(..., sep = "+"), 
                             as.list(model@response))
      return(pd)
    })
    plotData <- do.call(rbind, plotData)
    responseLevels <- sapply(model@models, function(m) m@response)
    plotData$response <- factor(plotData$response, levels = responseLevels)
    plotData$variable <- as.factor(plotData$variable)
    allLevels <- unique(c(levels(plotData$response), levels(plotData$variable)))
    variableLevels <- allLevels[allLevels %in% levels(plotData$variable)]
    plotData$variable <- factor(plotData$variable, levels = variableLevels)
    responseLevels <- allLevels[allLevels %in% levels(plotData$response)]
    plotData$response <- factor(plotData$response, levels = responseLevels)
    xLabel <- processData(model@models[[1]])@positionVar
    return(plotData)
  }
  .local(model, ...)
}

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
