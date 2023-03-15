## Analyses

This subdirectory contains the code to run the analyses from the paper. For Quarto (.qmd) files, the html rendered version is also available. The different parts can be run in any order, but they are numbered to follow the structure of the article.

-   `01_hawkes_process_simulation_example/` contains the code to simulate and plot a simple Hawkes process.
-   `02_simulation_interevent_times/` contains the code to simulate data in different conditions and infer interactions with an inter-event times method (Murphy et al., 2021). The results from the article were generated using `slurm` on a computing cluster (for this work, the computing facilities of the CC LBBE/PRABI was used). The code can also be run locally, but the full loop would take around 1400 hours to run.
    -   `01_launch_jobs.sh` loops `compute_inter_event_times.slurm` to generate the simulations in different conditions. To run this script locally, you should set the variable `run_locally` to `true` in this script and navigate to `analyses/02_simulation_interevent_times/`. Then, run `./launch_jobs.sh`. If run locally, the results will be written in the directory `outputs/02_simulation_interevent_times/simu_local/` (this name can be changed in `01_simu_infer.R`).
    -   `01-1_compute_inter_event_times.slurm` launches `01_simu_infer.R` with different parameters (locally or via `slurm`)
    -   `01-2_simu_infer.R` generates a single simulation and associated permutations
    -   `02_analyze_interevent_times.qmd` analyzes the results from the simulations
-   `03_simulation_hawkes/` contains the code to simulate data in different conditions and infer interactions with Hawkes processes.
-   `04_plot_map/` contains the code to plot the sites map.
-   `05_example_real_data/` contains the code to infer a Hawkes process from camera trap data.
    -   `01_prepare_data.R` pre-processes the data
    -   `02_infer.qmd` infers the Hawkes process
-   `06_circadian_rhythm/` contains the code to simulate camera trap data for species with circadian rhythms and infer a Hawkes process from those data (in Appendix S1).

## References

Murphy, A., Diefenbach, D. R., Ternent, M., Lovallo, M., & Miller, D. (2021). Threading the needle: How humans influence predator--prey spatiotemporal interactions in a multiple-predator system. Journal of Animal Ecology, 90(10), 2377--2390. <https://doi.org/10.1111/1365-2656.13548>
