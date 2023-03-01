## Code

This subdirectory contains the code to run the analyses from the paper.

For Quarto (.qmd) files, the html rendered version is also available.
The high-quality figures used in the article can be found in the `figures` sub-directory of each folder.

+ `circadian_rhythm` contains the code to simulate camera trap data for species with circadian rhythms and infer a Hawkes process from those data (in Appendix S1).
+ `example_real_data` contains the code to infer a Hawkes process from camera trap data.
  + `01_prepare_data.R` pre-processes the data (written in the `out` subfolder).
  + `02_infer.qmd` infers the Hawkes process.
+ `hawkes_process_simulation_example` contains the code to simulate and plot a simple Hawkes process.
+ `plot_map` contains the code to plot the sites map.
+ `simulation_hawkes` contains the code to simulate data in different conditions and infer interactions with Hawkes processes.
  + pre-computed results can be found in the `out` subfolder (the code takes about 1 hour to run).
+ `simulation_interevent_times` contains the code to simulate data in different conditions and infer interactions with an inter-event times method (Murphy et al., 2021).
  + `01_simu_infer.R` generates a single simulation and associated permutations.
  + `02_analyze_interevent_times.qmd` analyzes the results from the simulations.
  +  `compute_inter_event_times.slurm` launches `01_simu_infer.R` on a computer cluster using `slurm` (for this work, the computing facilities of the CC LBBE/PRABI was used).
  + `launch_jobs.sh` loops `compute_inter_event_times.slurm` to generate the simulations in different conditions.
  + The pre-computed results can be found in the `out` subfolder: 
    + `simu_interevents_cluster` contains all inference repetiitons with simulated data (that were run on the cluster) 
    + `out/pval_precomputed.csv` contains the pre-computed p-values associated to each interaction for each repetition (else, it takes a few minutes to compute).






## References

Murphy, A., Diefenbach, D. R., Ternent, M., Lovallo, M., & Miller, D. (2021). Threading the needle: How humans influence predator–prey spatiotemporal interactions in a multiple-predator system. Journal of Animal Ecology, 90(10), 2377–2390. https://doi.org/10.1111/1365-2656.13548
