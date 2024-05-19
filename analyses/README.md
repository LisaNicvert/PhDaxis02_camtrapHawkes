## Analyses

This subdirectory contains the code to run the analyses from my thesis. Subfolders 01 to 07 contain analyses for the linear multivariate Hawkes pricess (MHP) (chapter 3.2), and subfolder 08 contains the analyses for the non-linear MHP (chapter 3.3).

For Quarto (.qmd) files, the html rendered version is also available. The different parts can be run in any order, but they are numbered to follow the structure of the article.

-   `run_all.R` is a script to launch all analysis scripts for the article (PhD chapter 3.2) in a sequential order (except `02_simulation_interevent_times/01_launch_jobs.sh` because it is impractical on a personal computer; see below). **Warning: some analyses take time** (especially `03_simulation_hawkes/` that takes about 1 hour to run).
-   `01_hawkes_process_simulation_example/` contains the code to simulate and plot a simple Hawkes process.
-   `02_simulation_interevent_times/` contains the code to simulate data in different conditions and infer interactions with an inter-event times method (Murphy et al., 2021). This folder also contains some of the code run for Appendix S1: Section S1. The results from the article were generated using `slurm` on a computing cluster (for this work, the computing facilities of the CC LBBE/PRABI was used). To run the script on the LBBE/PRABI cluster, we used `Singularity` with the Docker provided in this repository. The code can also be run locally, but the full loop would take around 1400 hours to run.
    -   `01_launch_jobs.sh` simulates data in different conditions using `01-1_compute_inter_event_times.slurm`. To run `01_launch_jobs.sh` locally, you should set the variable `run_locally` to `true` in this script (**warning: if you run this analysis locally, use a smaller set of parameters or it will take too much time**). Navigate to `analyses/02_simulation_interevent_times/` and then run `./01_launch_jobs.sh`. If run locally, the results will be written in the directory `outputs/02_simulation_interevent_times/simu_local/` (this name can be changed in `01-2_simu_infer.R`).
    -   `01-1_compute_inter_event_times.slurm` launches `01-2_simu_infer.R` with different parameters (locally or via `slurm` with `Singularity`)
    -   `01-2_simu_infer.R` generates a single simulation and associated permutations
    -   `02_analyze_interevent_times.qmd` analyzes the results from the simulations.
-   `03_simulation_hawkes/` contains the code to simulate data in different conditions and infer interactions with Hawkes processes. This code takes about 1 hour to run when `perform_simu` and `perform_inference` are set to `TRUE` (as in the default setting). This folder also contains some of the code run for Appendix S1: Section S1.
-   `04_plot_map/` contains the code to plot the sites map.
-   `05_example_real_data/` contains the code to infer a Hawkes process from camera trap data. This folder also contains the code run for Appendix S1: Section S2 and the code to reproduce the figure from the introduction of part 3 in the thesis (capture event barplot by species).
    -   `01_prepare_data.R` pre-processes the data
    -   `02_infer.qmd` infers the Hawkes process
-   `06_circadian_rhythm/` contains the code to simulate camera trap data for species with circadian rhythms and infer a Hawkes process from those data (in Appendix S1: Section S3).
-   `07_bins_width/` contains the code to evaluate the importance of bins width on the inference (in Appendix S1: Section S4).
-   `08_ppstat/` contains the code to simulate, evaluate and infer datza with the non-linear MHP (chapter 3.3).

## References

Murphy, A., D. R. Diefenbach, M. Ternent, M. Lovallo, and D. Miller. 2021. Threading the needle: How humans influence predator–prey spatiotemporal interactions in a multiple-predator system. Journal of Animal Ecology **90**:2377–2390. <https://doi.org/10.1111/1365-2656.13548>
