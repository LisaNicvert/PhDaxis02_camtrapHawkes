## Outputs

This subdirectory contains the output of the corresponding analyses found in `analyses/`.

-   `02_simulation_interevent_times/` contains the outputs of the inference with an inter-event times method (Murphy et al., 2021).
    -   `simu_interevents_cluster/` contains all inference repetitions with simulated data (that were run on the computing facilities of the CC LBBE/PRABI)
    -   `out/pval_precomputed.csv` contains the pre-computed p-values associated to each interaction for each repetition (else, it takes a few minutes to compute).
-   `03_simulation_hawkes/` contains the output of the inference with Hawkes processes (the code takes about 1 hour to run).
    -   `inferenceBVL.rds` contains the inferred interactions
    -   `simulationsBVL.rds` contains the simulated datasets
-   `04_plot_map/` contains the shapefiles needed to plot the map (extracted from OpenStreetMap).
    -   `reserves/` contains the shapefiles for the reserves.
    -   `southern_africa/` contains the shapefiles for the countries borders.
-   `05_example_real_data/` contains the cleaned data for the inference of a Hawkes process from camera trap data.
