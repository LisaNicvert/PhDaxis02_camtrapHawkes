#!/bin/bash

varstrength=$1 # get first variable given as an argument to sbatch (strength)
vardtmax=$2 # get second variable (Tmax)
varrep=$3 # get third variable (repetition)
vari=$4 # get job index
varlocal=$5 # run locally?

begin=`date +%s`
currentDate=`date`

echo "=================================="
echo "DateTime: ${currentDate}"
echo "=================================="

echo "strength=${varstrength}"
echo "Tmax=${vardtmax}"
echo "rep=${varrep}"
echo "i=${vari} (seed)"
echo "run_locally=${varlocal}"


if ! $varlocal ; then
  # --- Run on the cluster
  cd /beegfs/home/lnicvert/camtrapHawkes
  /beegfs/data/soft/singularity3.11.0/bin/singularity run docker://gitlab-registry.in2p3.fr/stephane.delmotte/rocker_lbbe:master Rscript analyses/02_simulation_interevent_times/01-2_simu_infer.R $varstrength $vardtmax $varrep $vari $varlocal
else
  # --- Run locally 
  Rscript 01-2_simu_infer.R $varstrength $vardtmax $varrep $vari $varlocal
fi

end=`date +%s`
elapsed=$(echo "scale=2; ($end - $begin) / 60" | bc -l)

echo "Job finished"
echo "Time taken: ${elapsed} min"

