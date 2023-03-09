#!/bin/bash

# ========== Simulation parameters ==========
# --- Initialize strength list
strengthlist=(0.01 0.1 0.2 0.5 1)

# --- Initialize Tmax list
tmaxlist=(20 100 300 400 500)

# --- Initialize nrep (used for R seed)
nrep=30


# ========== sbatch parameters ==========
# --- Create directories (if they don't exist)
# Create error directory 
if [ ! -d out/sbatch/error ]; then
mkdir -p out/sbatch/error
fi
# Create output directory 
if [ ! -d out/sbatch/output ]; then
mkdir -p out/sbatch/output
fi

# --- Time
time="10:00:00"

# ========== Launch jobs ==========
i=734
# --- Iterate through strength
for strength in ${strengthlist[@]}
do
	# --- Iterate through Tmax
	for tmax in "${tmaxlist[@]}"
	do
		# --- Repeat permutation n times
		for rep in $(seq 15 1 $nrep)
		do
			i=$((i+1))
			sbatch --job-name str$strength-tmax$tmax-rep$rep-i$i --output out/sbatch/output/output_i$i-str$strength-tmax$tmax-rep$rep.out --error out/sbatch/error/error_i$i-str$strength-tmax$tmax-rep$rep.out --time=$time --nodes=1 --cpus-per-task=4 --mem=1G /beegfs/home/lnicvert/camtrap_hawkes/code/inter_event_times/compute_inter_event_times.slurm $strength $tmax $rep $i
		done 
	done
done

exit 0
