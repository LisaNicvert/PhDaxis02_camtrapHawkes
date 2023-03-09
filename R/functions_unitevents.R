# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-10-28
#
# Script Description: functions related to UnitEvents data handling, plotting, inference etc.

# Model creation ----------------------------------------------------------

#' Create background rates
#' 
#' Create a background rates vector.
#'
#' @param spont spontaneous coefficient 
#'   if unique value will be repeated
#'   else will be in the order provided
#' @param spp_names names of species
#'
#' @return Returns a matrix of dim (nspecies, 1) containing named spont rates
#' @export
create_spont_rates <- function(spont, spp_names){
  nspecies <- length(spp_names)
  
  # Spontaneous part
  if(length(spont) == 1){
    res <- rep(spont, nspecies)
  }else{
    if(length(spont) != length(spp_names)){
      print("Warning, the length of spont should be the same as the length of spp_names")
    }
    res <- as.numeric(spont)
  }
  
  names(res) <- spp_names
  
  res <- matrix(data = res,
                nrow = nspecies)
  rownames(res) <- spp_names

  return(res)
} 

#' Create interactions
#' 
#' Create an array of null interaction functions
#'
#' @param spp_names names of species
#' @param times times vector
#'
#' @return A (nspecies, nspecies, 1) array containing in each
#'   cell a (2, length(times)) matrix: top row is a vector of zero (null)
#'   interaction function), borrom row is the times vector.
#' @export
create_interactions <- function(spp_names, times){

  nspecies <- length(spp_names)

  # Create zero function
  funczero <- rep(0, length = length(times))
  
  funczero_unit <- matrix(c(times, funczero),
                          nrow = 2, byrow = TRUE)
  funczero_ue <- array(data = list(funczero_unit), 
                       dim = c(nspecies, nspecies,1)) 
  
  # Name rows and columns
  rownames(funczero_ue) <- spp_names
  colnames(funczero_ue) <- spp_names
  
  return(funczero_ue)
}

#' Create function shape
#' 
#' Create a function shape.
#'
#' @param funcshape exp, gamma or linear
#' @param times times vector
#' @param t half-life (time value for which the function reaches strength/2).
#'   only used when funcshape == "exp"
#' @param strength function max strength coefficient
#'
#' @return A numeric vector of the same length as the times vector containing 
#'   function values computed for each time.
#' @export
create_funcshape <- function(funcshape, times, t = 0.5, strength){

  if(funcshape == "exp"){
    # Create non-zero interaction function
    fs <- strength*exp(log(1/2)/t*times)
    
  }else if(funcshape == "gamma"){
    sh <- 7.5
    fs <- strength*stats::dgamma(times, shape = sh, scale = 1) # 1/(sh - 1)
    
  }else if(funcshape == "linear"){
    fs <- seq(strength, 0, length = length(times))
  }
  
  return(fs)
}


# True/false positives/negatives ------------------------------------------

#' Get bins sum
#' 
#' Get the sum of the bins for a Hawkes model.
#'
#' @param M an interaction matrix(K*K*Ncomp array)
#'   (K species and Ncomp different comportments (often Ncomp =1))
#'
#' @return A K*K matrix containing 1 if the interaction was inferred, else 0
get_bins_sum <- function(M) {
  
  res <-apply(M[,,1], c(1,2), FUN = function(m) sum(abs(m[[1]][2,])))
  res[res != 0] <- 1
  
  return(res)
}

#' Get bins sums
#'
#' Get the sum of the bins for a list of Hawkes model.
#' 
#' @param Mlist Mlist: a list of interaction matrix(K*K*Ncomp array)
#'   (K species and Ncomp different comportments (often Ncomp =1))
#'
#' @return A K*K matrix containing the number of times each value was inferred
get_bins_sum_list <- function(Mlist) {
  
  res <- matrix(data = 0,
                nrow = nspecies,
                ncol = nspecies)
  for (M in Mlist) {
    r <- get_bins_sum(M[[est]]$I)
    res <- res + r
  }
  return(res)
}

#' Compute true/false positives/negatives
#' 
#' Compute:
#'     TP (true positives)
#'     FP (false positives)
#'     TN (true negatives)
#'     FN (false negatives)
#'
#' @param Msimul a K*K*Ncomp array (K species Ncomp different comportments (often Ncomp =1))
#' @param Mtrue true array (K*K*Ncomp array)
#'
#' @return Returns a named vector with four elements named TP, FP, TN, FN
#' @export
compute_pos_neg <- function(Msimul, Mtrue){
  
  Msim_grouped <- get_bins_sum(Msimul)
  Mtrue_grouped <- get_bins_sum(Mtrue)
  
  # -- Positives
  P = Msim_grouped[Mtrue_grouped != 0]
  # True Positives
  TP = length(P[P != 0])
  
  # False negatives
  FN = length(P[P == 0])
  
  # -- negatives
  Neg = Msim_grouped[Mtrue_grouped == 0]
  
  # True negatives
  TN = length(Neg[Neg == 0])
  
  # False positives
  FP = length(Neg[Neg != 0]) 
  res = c(TP, FP, TN, FN)
  
  names(res) = c('TP', 'FP', 'TN', 'FN')
  return(res)
}

#' Compute several true/false positives
#' 
#' Evaluates several inferred models (in Msim_list) compared to the true
#' model used to generate data (Mtrue).
#' Computes:
#'     TP (true positives)
#'     FP (false positives)
#'     TN (true negatives)
#'     FN (false negatives)
#' @param Msim_list a list, each element has elements $BL and $BOL and
#'     in each there are $S (spontaneous part) and $I 
#'     (interaction coefficients, lists.)
#' @param Mtrue true model: array (K*K*Ncomp array)
#' @param est estimator to use (BL, BVL or BOL)
#'
#' @return Returns a df with columns TP, FP, TN, FN and as many rows as
#'   there are models in Msim_list.
#' @export
compute_pos_neg_list <- function(Msim_list, Mtrue, est = "BL"){
  
  N <- length(Msim_list)
  
  res <- c()
  for(i in 1:N){
    Msim <- Msim_list[[i]][[est]]$I
    
    pos_neg <- compute_pos_neg(Msimul = Msim, Mtrue = Mtrue)
    res <- rbind(res, pos_neg)
  }
  return(res)
}


# Get simulated models ----------------------------------------------------

#' Get models dataframe
#'
#' Get a list of all models in inferred_models plus the true model
#' 
#' @param inferred_models a list with 2 components
#'    $reinfer: list of models, each one inferred from a different dataset
#' which is assumed to have been generated by the same true model.
#'    $reinfer_parameters: metadata about the inferred models
#' @param models a list of true models to match the true model that have been 
#'   used to generate data from which models in $reinfer were inferred.
#' @param est estimator to use ($BL or $BOL)
#'
#' @return A dataframe  of models
#' @export
get_models_df <- function(inferred_models, 
                          models, 
                          est = "BL"){
  
  # Get the true model from which simulated data were inferred
  if("S" %in% names(models) & "I" %in% names(models) & length(models) == 2){ # Then it's only one model
    mtrue <- models
  }else if("modelname" %in% names(inferred_models$reinfer_parameters)){
    # Model name is directly available
    modelname <- inferred_models$reinfer_parameters$modelname
    mtrue <- models[[modelname]]
  }
  
  # Get the inferred models (without parameters)
  minf <- inferred_models$reinfer

  # Get the inferred models
  res <- list()
  for(i in 1:(length(minf))){
    res[[i]] <- minf[[i]][[est]]
  }
  
  res[[length(res) + 1]] <- mtrue
  
  # Name the models list (s + repetition #, and "true" for the true model)
  names(res) <- c(paste0("s", seq(1,length(res)-1)),
                  "true")
  
  # Transform list to dataframe
  res <- mlist_to_df(res)
    
  return(res)
}

#' List of models to dataframes
#' 
#' Transform a list of UnitEvent models to a dataframe
#' 
#' @param mlist list of models of type UnitEvents 
#'   each component of mlist must be a list with 2 named 
#'   elements $S and $I (corresponds to a model$BL, BOL or $BVL output of BoxLasso)
#' @param comp the comportment number to extract
#'
#' @return A dataframe with 4 or 5 columns:
#'   time: time for the interaction functions
#'   excitefunc: value of the excitation function
#'   from: from species
#'   to: to species
#'   spont: spontaneous rate of the "to" species
#'   rep (if mlist is a named list): name of the i-th element of mlist
#' @export
mlist_to_df <- function(mlist, comp = 1){
  
  mli <- mlist
  
  # Name list if not named already
  if(is.null(names(mli))){
    names(mli) <- paste0("L", 1:length(mli))
  }
  
  mall <- lapply(mli, ue_model_to_df, comp = comp)
  
  # Add a column  with repetition ID
  mall <- lapply(seq_along(mall), 
                 function(i){
                   mall[[i]] %>% 
                     mutate(rep = rep(names(mall)[i],
                                      nrow(mall[[i]])))})
  
  names(mall) = names(mli) # Add names
  
  # Get one bigass df
  res <- do.call("rbind", mall)
  
  return(res)
}

# Compute rates -------------------------

#' Compute rates
#' 
#' Compute the rate corresponding to data with a given model.
#' 
#' @param model a dataframe corresponding to the model estimated with UnitEvents.
#'   This dataframe can be obtained with the function "ue_model_to_df".
#'   Must have columns time, excitefunc, to, from and spont.
#' @param data occurrence data. Must have one column stamp and one column species. 
#' Currently only one camera is supported.
#' @param timestep timestep for the function discretisation (days). 
#' Should be smaller than delta used in the model for good results. Defaults to 0.01.
#'
#' @return A dataframe with columns time, lambda and species.
#'   time: the time.
#'   lambda: the intensity function at time t.
#'   species: the species for which the intensity is computed.
#' @export
compute_rate <- function(model, data, timestep = 0.01){
  
  # Occurrence times
  stamps <- data$stamp
  
  # --- Times vector
  # goes up to (including) timestep + function support
  time <- seq(min(stamps), max(stamps) + max(model$time) + timestep, 
              by = timestep)
  
  # --- Compute the rate for each species
  species_list <- unique(model$from)
  
  # Initialise df
  allres <- data.frame()
  
  for(target in species_list){ # for each species
    
    # Spont rate for the target species
    spont <- unique(model$spont[model$to == target])
    
    # Initialise intensity
    lambda <- rep(spont, length(time))
    
    # Species
    species <- rep(target, length(time))
    
    # Dataframe for the target species
    res <- data.frame(time, lambda, species)
    
    for(occ in 1:nrow(data)){ # Get all occurrences of all species
      species <- data$species[occ]
      t <- data$stamp[occ]
      
      # Get the function species -> target
      func <- model %>% 
        filter(to == target & from == species) %>%
        select(time, excitefunc) %>%
        arrange(time)
      funcval <- func$excitefunc
      
      for(k in 1:(nrow(func)-1)){
        # get k-th bin
        fk <- func$excitefunc[k]
        # Get time span on which k-th bin acts
        tmin <- t + func$time[k]
        tmax <- t + func$time[k+1]
        # Compute rate value between t and t + delta (delta = bin width)
        res$lambda[res$time > tmin & res$time <= tmax] <- res$lambda[res$time > tmin & res$time <= tmax] + fk
      }
    }
    
    # Bind rows for all target species
    allres <- allres %>% bind_rows(res)
    
  }
  
  return(allres)
  
}

# Inference parameters --------------------------

#' Compute box
#' 
#' Compute time windows covering the whole dataset for a given dataset d.
#'
#' @param d dataframe with columns stamp, datetime and cameraID
#' @param use_stamps use timestamps or dates?
#'
#' @return Returns the extimated box (matrix)
#'   line 1: Camera number
#'   line 2: start time
#'   line 3: stop time
#'   line 4: comportment ID
#'   each column corresponds to one camera.
#' @export
compute_box <- function(d, use_stamps = FALSE){
  
  # --- Compute stamps from custom origin
  if (!use_stamps) {
    orig <- min(d$datetime)
    hour(orig) <- 00
    minute(orig) <- 00
    second(orig) <- 00
    
    dstamp <- add_stamps(d, origin = orig)
  } else {
    dstamp <- d
  }
  
  Ntrial <- as.numeric(length(unique(d$cameraID)))
  # --- Summarise data
  Tmax.df <- d %>% group_by(cameraID) %>% 
    summarise(Tmin = min(stamp),
              Tmax = max(stamp)) %>%
    arrange(cameraID) # because in the matrix is in that order
  
  # Correct boxed to include bounds
  Tmax <- Tmax.df$Tmax + 0.01
  Tmin <- Tmax.df$Tmin - 0.01
  
  BoxEst <- matrix(c(1:Ntrial,
                     Tmin, 
                     Tmax, 
                     rep(1, Ntrial)),
                   nrow = 4, byrow=TRUE)
  # Name BoxEst columns according to the camera ID
  # (Useless, for better code readability)
  colnames(BoxEst) <- Tmax.df$cameraID
  
  return(BoxEst)
}

#' Get IDs
#' 
#' Get a unique row ID for each combination of species/camera.
#' 
#' @param df a dataframe that must have columns
#'     species
#'     cameraID
#' @param delete_missing_species delete species not seen at all cameras?
#'
#' @return a dataframe with columns
#'    cameraID
#'    species
#'    rowid (species_cameraID)
#' The df is arranged with all species for the same site first
#' 
#' @export
get_ids <- function(df, delete_missing_species = FALSE){
  
  if(delete_missing_species){ # if asked first don't give an ID for species
    # which are not seen at all sites
    sites_count <- length(unique(df$cameraID))
    filtered <- df %>% group_by(species) %>% 
      summarise(n_sites = length(unique(cameraID)))
    spp <- filtered %>% filter(n_sites == sites_count) # get only species present on all 4 cameras
    df_filtered <- df %>% filter(species %in% spp$species)
  }
  
  # Get species and site list
  species_list <- unique(as.character(sort(df$species)))
  sites_list <- unique(as.character(sort(df$cameraID)))
  
  # Initialise empty vectors
  ids_df <- expand.grid(sites_list, species_list)
  colnames(ids_df) <- c("cameraID", "species")
  
  # Add ID
  ids_df <- ids_df %>%
    mutate(rowid = paste(species, cameraID, sep = "_")) %>%
    select(rowid, species, cameraID)
  
  # Rearrange
  ids_df <- ids_df %>% arrange(species, cameraID)
  
  return(ids_df)
}

# Inference --------------------------

#' Infer a Hawkes model
#' 
#' Infer Hawkes model for data d.
#'
#' @param d dataset to infer from. Must have columns:
#'   snapshotName
#'   cameraID
#'   stamp
#'   datetime
#'   count
#' @param k number of bins for the inference function
#' @param delta binwidth
#' @param Z cutoff
#' @param gamma penalization LASSO parameter
#' @param scale scale for precision of calculations (defaults to 10000)
#' @param use_stamps use timestamps or dates?
#'
#' @return Returns the output in the same form as BoxLasso() function of UnitEvents
#' A list with 2 (if n(species) = 1) or 3 elements (if n(species) > 1)
#'   $BL
#'   $BVL
#'   $BOL
#' @export
infer <- function(d, 
                  k = 12, 
                  delta = 2/24, 
                  Z = 0,
                  gamma = 0.5, 
                  scale = 10000,
                  use_stamps = FALSE){

  # --- Inference parameters from data
  Nneur = as.numeric(length(unique(d$snapshotName)))
  Ntrial = as.numeric(length(unique(d$cameraID)))
  
  BoxEst <- compute_box(d = d,
                        use_stamps = use_stamps)
  
  # Rename snapshotName -> species to fit for more generic function df_to_matrix
  d <- d %>% rename("species" = "snapshotName")
  
  # Then transform data to suitable form for analysis
  M <- df_to_matrix(d)
  
  # Lasso inference
  inf <- BoxLasso(DataNeur = M,
                  Nneur = Nneur,
                  scale = scale,
                  del = delta, 
                  k = k,
                  Z = Z,
                  BoxEst = BoxEst,
                  gamma = gamma)
  
  return(inf)
}

#' Format the inference output
#' 
#' Transform output of inference into a list of coefficient matrices.
#' 
#' @param inf the output of the infer or BoxLasso. Must be a list with 3 elements:
#'   $BL
#'   $BVL
#'   $BOL
#' @param species_names a list of species names to rename rows and columns (optional)
#' @param k number of bins
#' @param delta binwidth
#'
#' @return A list of length 3:
#'   $BL
#'   $BVL
#'   $BOL
#' Each element is a list of length 2 (and is the named output of the UnitEvent's
#'   function coeff2interac):
#'   $S (matrix (n(species), Ncomp)). The matrix's rows are named
#'     with species_names.
#'   $I (array (n(species), n(species), Ncomp)). The array's first
#'   and second dimensions are named with species_names.
#' @export
unpack_inf <- function(inf, species_names = NA, k, delta){
  
  bench = list() # initialise list
  bench[[1]] = coeff2interac(inf$BL, k, delta) # lasso estimators
  bench[[2]] = coeff2interac(inf$BOL, k, delta) # least squares estimators with LASSO support
  bench[[3]] = coeff2interac(inf$BVL, k, delta) # estimator for which some coefficients vanish
  names(bench) <- c("BL", "BOL", "BVL")
  
  if(length(species_names) == 1){ # it's maybe a default NA
    if(is.na(species_names)){ # it's the default
      species_names <- seq(1, nrow(bench[[1]]$S)) # Replace species names with 1..N
    }
  }else{
    species_names <- sort(species_names)
  }
  
  # Rename rows and columns
  for (i in 1:length(bench)){
    row.names(bench[[i]]$S) <- species_names
    
    row.names(bench[[i]]$I) <- species_names
    colnames(bench[[i]]$I) <- species_names
  }
  
  return(bench)
}

#' Reinfer Hawkes model
#' 
#' Reinfer Hawkes model from a dataset generated with a Hawkes model.
#' 
#' @param M dataset to infer from, a DataNeur matrix (output of HawkesMulti
#'   from the UnitEvents package). It is an array (nspecies, ncameras, nmax+1)
#'   where nmax is the maximum number of occurrences for one species on
#'   one camera.
#' @param Ntrial number of cameras
#' @param k number of bins for the interaction functions
#' @param delta binwidth
#' @param gamma penalization LASSO parameter
#' @param scale scale for precision of calculations (defaults to 10000)
#'
#' @return Returns the output in the same form as BoxLasso() function of UnitEvents
#' A list with 2 (if n(species) = 1) or 3 elements (if n(species) > 1)
#'   $BL
#'   $BVL
#'   $BOL
#' @export
reinfer <- function(M, 
                    Ntrial,
                    k = 12, 
                    delta = 2/24, 
                    gamma = 0.5, 
                    scale = 10000){
  
  # Compute Tmax for each camera for the boxes later
  d <- ue_to_df(M)
  Tmax.df <- d %>% group_by(camera) %>% 
    summarise(Tmax = max(stamp)) %>%
    arrange(camera) # because in the matrix is in that order
  Tmax <- Tmax.df$Tmax
  
  # Inference parameters from data
  Nneur = nrow(M)/Ntrial
  
  if(Ntrial > 1){
    BoxEst <- matrix(c(1:Ntrial,
                       rep(0,Ntrial), 
                       Tmax, 
                       rep(1, Ntrial)),
                     nrow = 4, byrow=TRUE)
  }else{
    BoxEst <- matrix(c(1,
                       0, 
                       Tmax, 
                       1),
                     nrow = 4, byrow=TRUE)
  }
  
  # Lasso inference
  inf <- BoxLasso(DataNeur = M,
                  Nneur = Nneur,
                  scale = scale,
                  del = delta, 
                  k = k,
                  BoxEst = BoxEst,
                  gamma = gamma)
  
  return(inf)
}

# Convert data formats --------------------------

#' Transform dataframe to matrix
#' 
#' Transforms the dataframe into a matrix ready for Lasso inference.
#' 
#' @param df A dataframe with columns
#'   cameraID
#'   stamp
#'   species
#'   count
#'
#' @return an array (nspecies, ncameras, nmax+1) in the same format as the output 
#'   of HawkesMulti from the UnitEvents package). nmax is the maximum number 
#'   of occurrences for one species on one camera.
#' @export
df_to_matrix <- function(df){
  
  ids_df <- get_ids(df, delete_missing_species = FALSE)
  
  # Remove useless columns 
  prepmat <- df %>% select(stamp, species, cameraID) %>% 
    mutate(species = as.character(species)) %>%
    mutate(cameraID = as.character(cameraID))
  
  # Create an ID for species and camera
  prepmat <- prepmat %>% arrange(species, cameraID, stamp) %>% # arrange by species then camera  
    # then by timestamp (because later we want data grouped by species then camera)
    mutate(rowid = paste(species, cameraID, sep = "_"))
  
  # Transform table: each set of times is stored as a list and its length
  # is stored in col len.
  prepmat <- prepmat %>% ungroup() %>% group_by(rowid) %>% 
    summarise(times = list(sort(stamp)), 
              len = n())
  
  # Add missing row IDs (in case one species was not seen at one camera, because
  # we need the same trial numbeer for each species)
  prepmat_full <- ids_df %>% left_join(prepmat, by = 'rowid') 
  
  # Re-sort per species, then camera ID (should be useless but just in case)
  prepmat_full <- prepmat_full %>% arrange(species, cameraID) 
  
  maxlen <- max(prepmat_full$len, na.rm = TRUE) # get future column number
  
  prepmat_full <- prepmat_full %>% mutate(rep_zero = maxlen - len) # get the number of zeroes to be repeated
  
  M <- matrix(nrow = nrow(prepmat_full), ncol = maxlen)
  
  for(i in 1:nrow(prepmat_full)){
    if(!is.na(prepmat_full$rep_zero[i])){ # case there really were observations
      if(prepmat_full$rep_zero[i] != 0){
        M[i,] <- c(prepmat_full$times[[i]], rep(0, prepmat_full$rep_zero[i]))
      }else{
        M[i,] <- prepmat_full$times[[i]]
      }
    }else{ # no observations => fill with zeros
      M[i,] <- rep(0, maxlen) # fill the row with zeros
    }
  }
  
  # Rename rows
  row.names(M) <- prepmat_full$rowid
  
  # Add first column
  M <- counting_spikes(M)
  colnames(M) <- c("nspikes", 1:(ncol(M)-1)) # rename first column
  
  return(M)
}

#' Model to dataframe
#' 
#' Transforms the specification of a UnitEvents model 
#' to a dataframe.
#'
#' @param mod  a list of length 2 (and is the named output of the UnitEvent's
#'   function coeff2interac):
#'   $S (matrix (n(species), Ncomp)).
#'   $I (array (n(species), n(species), Ncomp)).
#' @param comp the comportment to extract
#'
#' @return a dataframe with columns
#'   time: the time (support of interaction functions)
#'   excitefunc: values of the functions for a corresponding time t
#'   from: species "from" which the interaction is
#'   to: species "towards" which the interaction is
#' @export
ue_model_to_df <- function(mod, comp = 1){

  # Get the interaction part
  mod_interacfunc <- mod$I[,, comp]
  
  # Get spontaneous part 
  mod_spont <- as.matrix(mod$S[, comp], ncol = 1)
  
  # --- Interaction
  # Transform list to dataframe
  mdf <- lapply(mod_interacfunc, as.data.frame)
  
  # Get index corresponding to each species
  nspp <- nrow(mod_spont)
  species_names <- rownames(mod$S)
  species <- cbind(rep(species_names, each = nspp), 
                   rep(species_names, nspp))
  colnames(species) <- c("from", "to")
  
  # transpose dataframes to get 1st column = times and 2nd column = funcion value
  mdf <- lapply(mdf, t)
  
  # Add species id
  nstep <- nrow(mdf[[1]]) # number of steps (repeats for species)
  mdf.spp <- lapply(seq_along(mdf), 
                    function(i){cbind(mdf[[i]], 
                                      matrix(rep(species[i,], nstep), 
                                             ncol = 2, 
                                             byrow = TRUE))})
  
  mdf.all <- do.call(rbind, mdf.spp)
  
  colnames(mdf.all) <- c("time", "excitefunc", "from", "to")
  
  # Add spont part
  mod_spont <- mod_spont %>% as.data.frame %>% 
    tibble::rownames_to_column("to")
  
  # Merge
  res <- data.frame(mdf.all) %>% left_join(mod_spont, by = 'to') %>%
    rename("spont" = "V1")
  # Convert to numeric
  res$spont <- as.numeric(res$spont)
  res$excitefunc <- as.numeric(res$excitefunc)
  res$time <- as.numeric(res$time)
  return(res)
}


#' Model to graph
#' 
#' Adaptation of plot_graph_Hawkes function in UnitEvents, 
#' that returns the igraph object instead of plotting it.
#'
#' @param IS a list of length 2 (the output of the UnitEvent's function coeff2interac):
#'   $S (matrix (n(species), Ncomp)).
#'   $I (array (n(species), n(species), Ncomp)).
#' @param neurnames species_names
#'
#' @return a tbl_graph object representing the model.
#'   nodes: 1 column: names (with species names).
#'   edge: 3 columns: from, to, weight.
#' @export
ue_to_graph <- function(IS, neurnames){
  
  spontComp = IS[[1]]
  interacComp = IS[[2]]
  
  Nneur = nrow(spontComp)
  Ncomp = ncol(spontComp)
  
  comp = 1
  mescolors = grDevices::colorRampPalette(c("red", "magenta", "gold", 
                                                 "green", "blue"))(Nneur)
  
  nul = which(spontComp[, comp] < 10^(-7))
  color = rep("red", Nneur)
  color[nul] = rep("cyan", length(nul))
  A = matrix(data = 0, ncol = Nneur, nrow = Nneur)
  B = A
  C = A
  for (i in 1:Nneur) {
    for (j in 1:Nneur) {
      part = interacComp[i, j, comp][[1]][1, ]
      valpart = interacComp[i, j, comp][[1]][2, ]
      lv = length(valpart)
      l1norm = sum(abs(valpart[1:(lv - 1)]) * diff(part))
      if (l1norm > 10^(-10)) {
        A[j, i] = l1norm
        B[j, i] = j
        C[j, i] = 1
      }
    }
  }
  for (i in 1:Nneur) {
    for (j in 1:Nneur) {
      if ((i != j) & (C[i, j] == 1) & (C[j, i] == 1)) {
        C[j, i] = 0.2
        C[i, j] = 0.2
      }
    }
  }
  g = igraph::graph.adjacency(A, weighted = TRUE, mode = "directed")
  if (Ncomp > 1) {
    for (comp in 2:Ncomp) {
      nul = which(spontComp[, comp] < 10^(-7))
      color = rep("red", Nneur)
      color[nul] = rep("cyan", length(nul))
      A = matrix(data = 0, ncol = Nneur, nrow = Nneur)
      B = A
      C = A
      for (i in 1:Nneur) {
        for (j in 1:Nneur) {
          part = interacComp[i, j, comp][[1]][1, ]
          valpart = interacComp[i, j, comp][[1]][2, 
          ]
          lv = length(valpart)
          l1norm = sum(abs(valpart[1:(lv - 1)]) * diff(part))
          if (l1norm > 10^(-10)) {
            A[j, i] = l1norm
            B[j, i] = j
            C[j, i] = 1
          }
        }
      }
      for (i in 1:Nneur) {
        for (j in 1:Nneur) {
          if ((i != j) & (C[i, j] == 1) & (C[j, i] == 
                                           1)) {
            C[j, i] = 0.2
            C[i, j] = 0.2
          }
        }
      }
      g = igraph::graph.adjacency(A, weighted = TRUE, mode = "directed")
    }
  }
  gr <- as_tbl_graph(g) %>% activate(nodes) %>% mutate(names = neurnames)
  return(gr)
}

#' Name data
#'
#' Attribute colnames and rownames to dataneur generated with Hawkesmulti.
#' 
#' @param dn an array (nspecies, ncameras, nmax+1) in the same format as the output 
#'   of HawkesMulti from the UnitEvents package). 
#' @param species_names ordered vector of species names.
#' @param cameras_names ordered vector of cameras names.
#'
#' @return The same matrix as the input (dn) but with named rows and columns.
#' @export
name_dataneur <- function(dn, species_names, cameras_names = NA){
  
  # copy data
  res <- dn
  
  # rename columns
  if(ncol(res) ==  1){ # case no observations
    colnames(res) <- "nspikes"
  }else{ # at least one observation
    colnames(res) <- c("nspikes",1:(ncol(dn)-1))
  }
  
  
  # get Ntrial
  Ntrial <- nrow(dn)/length(species_names)
  camnames <- paste0("C", 1:Ntrial)
  
  # rename rows
  nam <- sapply(species_names, FUN = function(s){
    r <- rep(s, Ntrial)
    paste(r, camnames, sep = "_")
  })
  names <- as.vector(matrix(nam, nrow = Ntrial*length(species_names)))
  row.names(res) <- names
  
  return(res)
}

#' Simolation to dataframe
#' 
#' Convert output of simulation generated with HawkesMulti to a dataframe 
#' (rows must be named).
#' 
#' @param ue an array (nspecies, ncameras, nmax+1) in the same format as the output 
#'   of HawkesMulti from the UnitEvents package). 
#'
#' @return Returns a dataframe with columns:
#'   species
#'   camera
#'   stamp
#' @export
ue_to_df <- function(ue){
  
  ID <- rownames(ue)
  IDs <- strsplit(ID, "_", fixed = TRUE)
  species <- sapply(IDs, function(i) i[1])
  cams <- sapply(IDs, function(i) i[2])
  
  simudf <- as.data.frame(ue[,2:ncol(ue)])
  
  simudf$species <- species
  simudf$camera <- cams  
  
  simudf <- simudf %>% select(species, camera, everything()) %>%
    tidyr::pivot_longer(cols = 3:ncol(simudf), values_to = "stamp")%>%
    select(-name)
  
  if(nrow(simudf) == 0){
    message("Warning: empty dataframe")
  }
  return(simudf)
}
