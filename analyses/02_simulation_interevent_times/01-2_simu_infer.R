# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-12-02
#
# Script Description: script to simulate datasets under different
# models, with different lengths, to infer interactions with
# inter-event-times.


# Libraries ---------------------------------------------------------------
# Install dependencies automatically (useful not to ask when script run with Rscript)
devtools::install_deps(upgrade = "never")
# Load functions from main folder
devtools::load_all()

# Load additional packages
packages <- c("here", "foreach", "doParallel", "UnitEvents",
              "magrittr", "dplyr", "tidyr")
base::lapply(packages, require)

registerDoParallel(cores = 4)
print(paste("getDoParWorkers:", getDoParWorkers()))

# Parameters --------------------------------------------------------------
# --- Model params
# strength_list <- c(0.01, 0.1, 0.2, 0.5, 1) # interaction strength (in args)
spont <- 0.1
nspecies <- 5
t <- 0.5 # half life
duration_days <- 2

# --- Simu params
# Tmax_list <- c(20, 100, 300, 400, 500) # Tmax (in args)
# nrep_each <- 30 # Number of datasets generated (in args)
n_cameras <- 30
perm_count <- 999 # Number of permutation

# --- Results writing
out_dir_cluster <- here("outputs/02_simulation_interevent_times/simu_interevents_cluster_test")
out_dir_local <- here("outputs/02_simulation_interevent_times/simu_local")

# Parse parameters --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# args <- c(0.1, 400, 1, 1, "true") # For test purposes if the script is not launched with launch_jobs.sh

# test if there is at least one argument: if not, return an error
if (length(args) != 5) {
  stop("You must supply exactly 5 arguments: strength, Tmax, rep, i and run_locally.")
}

strength <- as.numeric(args[1])
Tmax <- as.numeric(args[2])
rep <- as.numeric(args[3])
i <- as.numeric(args[4]) # used only for seed
local <- as.character(args[5])

# Set seed
set.seed(i)

# Define results folder name
if (local == "true") {
  out_dir <- out_dir_local
} else {
  out_dir <- out_dir_cluster
}

# Create results folder ---------------------------------------------------
if(!dir.exists(out_dir)) {
  print(paste("Create out_dir", out_dir))
  dir.create(out_dir, recursive = TRUE)
}

# Code --------------------------------------------------------------------

## Parameters for inter-event times ----------------------------------------
# Create variable response names
spp_names <- paste0("s", 1:nspecies)

# Get all edges to infer
all_pairs <- expand_grid(spp_names, spp_names)
colnames(all_pairs) <- c("from", "to")
all_pairs <- all_pairs %>% filter(from != to)

## Prepare true model ------------------------------------------------------
print("Prepare model ---------------------------------------")
# Times vector
times <- seq(0, duration_days, length.out = 200)

interac_null <- create_interactions(spp_names = spp_names,
                                    times = times)

# Spont rates
spont_rates <- create_spont_rates(spont, spp_names)

# Prepare function 
fs <- create_funcshape("exp", times = times, 
                       strength = strength)
fs[length(fs)] <- 0

func <- matrix(data = c(times, fs),
               nrow = 2, byrow = TRUE)

# Modify nullmodel where needed
interac <- interac_null
interac[2,1,1][[1]] <- func
interac[3,2,1][[1]] <- func
interac[4,2,1][[1]] <- func

# Store model in UE format
model <- list("S" = spont_rates, 
              "I" = interac)

## Simulate data -----------------------------------------------------------
print("Simulate data ---------------------------------------")
# Create BoxComp (only one comportment)
BoxComp <- array(c(list()), c(n_cameras,1))

for (j in 1:n_cameras){
  BoxComp[[j]] = matrix(c(0, Tmax, 1), nrow=3)
}

dat_simul <- Hawkesmulti(Ntrial = n_cameras, 
                         delay = 0, 
                         BoxComp = BoxComp, 
                         spontComp = model$S, 
                         interacComp = model$I)

dat_simul <- name_dataneur(dat_simul,
                           species_names = spp_names)

## Infer with inter-event times --------------------------------------------
print("Infer with inter-event times ---------------------------------------")

dat_simul_df <- ue_to_df(dat_simul) %>% 
  filter(stamp != 0)

comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}
  
res_list <- foreach(p = 1:nrow(all_pairs), 
                    .combine = "comb")  %dopar% {
  pair <- all_pairs[p,]
  # Separate df for each species
  sp_from <- dat_simul_df %>% filter(species == pair$from)
  sp_to <- dat_simul_df %>% filter(species == pair$to)
  
  ### Observed -------------
  # Compute observed interval
  intervals <- compute_intervals(sp_from = sp_from, sp_to = sp_to)
  
  # Compute median of min inter-event time
  intervals_median <- median(intervals, na.rm = TRUE)
  
  ### Permuted -------------
  # Add begin and end
  sp_from <- sp_from %>%
    group_by(camera) %>%
    mutate(Begin = min(stamp),
           End = max(stamp)) %>%
    ungroup()
  
  sp_to <- sp_to %>%
    group_by(camera) %>%
    mutate(Begin = min(stamp),
           End = max(stamp))
  
  # Compute permuted intervals
  permuted_values <- numeric(perm_count)
  for(s in 1:perm_count) {
    permuted_intervals <- compute_intervals_permute(sp_from = sp_from,
                                                    sp_to = sp_to)
    
    # Compute randomized median of min inter-event time
    permuted_intervals_median <- median(permuted_intervals, na.rm = TRUE)
    
    permuted_values[s] <- permuted_intervals_median
  }
    
  ### Store results -------------
  res <- data.frame("from" = pair$from, 
                    "to" = pair$to, 
                    "interval" = intervals_median, 
                    "rep" = rep,
                    "Tmax" = Tmax,
                    "strength" = strength)
  
  res_sim <- data.frame("from" = rep(pair$from, perm_count), 
                        "to" = rep(pair$to, perm_count), 
                        "interval" = permuted_values, 
                        "rep" = rep(rep, perm_count),
                        "Tmax" = rep(Tmax, perm_count),
                        "strength" = rep(strength, perm_count))
  
  res_list <- list("observed" = res,
                   "permuted" = res_sim)
}


# Write data --------------------------------------------------------------
print("Write data ---------------------------------------")

suffix <- paste(Tmax, 
                gsub(strength, 
                     pattern = ".", replacement = "", 
                     fixed = TRUE),
                rep, 
                i,
                sep = "_")

simname <- paste0("sim_", suffix, ".csv")
print(paste("Write sim to", file.path(out_dir, simname)))
write.csv(res_list$permuted, file.path(out_dir, simname),
          row.names = FALSE)

valname <- paste0("val_", suffix, ".csv")
print(paste("Write val to", file.path(out_dir, valname)))
write.csv(res_list$observed, file.path(out_dir, valname),
          row.names = FALSE)
