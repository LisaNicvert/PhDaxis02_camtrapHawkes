#!/bin/bash

varstrength=$1 # get first variable given as an argument to sbatch (strength)
vardtmax=$2 # get second variable (Tmax)
varrep=$3 # get third variable (repetition)
vari=$4 # get job index

begin=`date +%s`
currentDate=`date`

echo "=================================="
echo "DateTime: ${currentDate}"
echo "=================================="

echo "strength=${varstrength}"
echo "Tmax=${vardtmax}"
echo "rep=${varrep}"
echo "i=${vari} (seed)"

cd /beegfs/home/lnicvert/camtrap_hawkes

singularity run docker://gitlab-registry.in2p3.fr/stephane.delmotte/rocker_lbbe:master Rscript code/inter_event_times/01_simu_infer.R $varstrength $vardtmax $varrep $vari

end=`date +%s`
elapsed=$(echo "scale=2; ($end - $begin) / 60" | bc -l)

echo "Job finished"
echo "Time taken: ${elapsed} min"

