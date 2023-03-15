#!/bin/bash

# ========== Simulation parameters ==========
# --- Initialize strength list
strengthlist=(0.01 0.1 0.2 0.5 1)

# --- Initialize Tmax list
tmaxlist=(20 100 300 400 500)

# --- Initialize nrep (used for R seed)
nrep=30

# --- Run the script locally or on a computing cluster?
run_locally=false

if ! $run_locally ; then
  # ========== navigate to directory ==========
  cd /beegfs/home/lnicvert/camtrapHawkes/
  
  # ========== sbatch parameters ==========
  # --- Create directories (if they don't exist)
  # Create error directory 
  if [ ! -d outputs/02_simulation_interevent_times/sbatch/error ]; then
  mkdir -p outputs/02_simulation_interevent_times/sbatch/error
  fi
  # Create output directory 
  if [ ! -d outputs/02_simulation_interevent_times/sbatch/output ]; then
  mkdir -p outputs/02_simulation_interevent_times/sbatch/output
  fi
  
  # --- Time
  time="10:00:00"
fi

# ========== Launch jobs ==========
i=0
# --- Iterate through strength
for strength in ${strengthlist[@]}
do
	# --- Iterate through Tmax
	for tmax in "${tmaxlist[@]}"
	do
		# --- Repeat permutation n times
		for rep in $(seq 1 1 $nrep)
		do
			i=$((i+1))
			
			if ! $run_locally ; then
			  # Used on the cluster
			  sbatch --job-name str$strength-tmax$tmax-rep$rep-i$i --output outputs/02_simulation_interevent_times/sbatch/output/output_i$i-str$strength-tmax$tmax-rep$rep.out --error outputs/02_simulation_interevent_times/sbatch/error/error_i$i-str$strength-tmax$tmax-rep$rep.out --time=$time --nodes=1 --cpus-per-task=4 --mem=1G /beegfs/home/lnicvert/camtrapHawkes/analyses/02_simulation_interevent_times/01-1_compute_inter_event_times.slurm $strength $tmax $rep $i $run_locally
			else
			  # Used locally
			  ./01-1_compute_inter_event_times.slurm $strength $tmax $rep $i $run_locally
			fi
		done 
	done
done

exit 0
