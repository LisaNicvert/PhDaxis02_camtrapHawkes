# Analyze camera trap data with the Hawkes process

This repository contains the code to reproduce the analyses and figures from the article *Using the Hawkes model to study spatio-temporal interactions between multiple species from camera trap data*.

## Packages installation

In order to run the analyses, the R package `UnitEvents` must first be installed separately.

### Install `UnitEvents`
To install `UnitEvents` (Lambert et al., 2018), see instructions to install from source [here](https://sourcesup.renater.fr/frs/?group_id=3267).

### Dependencies

In order to automatically install dependencies, the R package `devtools` must be intalled. It can be installed with:

```r
install.packages("devtools")
```

Then, other needed packages will automatically be loaded at the beginning of each analysis script. 

The functions contained in `R/` are loaded with `devtools::load_all()`. 
Other needed packages are loaded (and if needed, installed) for each analysis script with a custom `camtrapHawkes::require` function.

## Contents

+ `analyses/` contains all code to run the analyses
+ `data/` contains the camera trap data and species silhouettes used for plotting
+ `figures/` contains the high-quality figures used in the article generated by the analyses
+ `man/` contains the documentation fir the functions written in `R/`
+ `outputs/` contains the outputs generated during the analyses
+ `R/` contains the R functions developed for this project
+ `DESCRIPTION` contains the project metadata (author, date, dependencies, etc.)
+ `NAMESPACE` contains the namespace information for the functions in the `R/` folder.

## References
Lambert, R. C., Tuleau-Malot, C., Bessaih, T., Rivoirard, V., Bouret, Y., Leresche, N., & Reynaud-Bouret, P. (2018). Reconstructing the functional connectivity of multiple spike trains using Hawkes models. Journal of Neuroscience Methods, 297, 9–21. https://doi.org/10.1016/j.jneumeth.2017.12.026

## Session info

Below is the output of the `sessionInfo()` for the machine that was used to run the code:

```r
> sessionInfo()
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8       
 [4] LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
 [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.4.1            lubridate_1.9.2          dplyr_1.1.0             
[4] here_1.0.1               camtrapHawkes_0.0.0.9000

loaded via a namespace (and not attached):
 [1] viridis_0.6.1       pkgload_1.2.4       tidyr_1.2.0         tidygraph_1.2.0    
 [5] viridisLite_0.4.0   splines_4.1.3       foreach_1.5.2       ggraph_2.1.0       
 [9] prodlim_2019.11.13  brio_1.1.3          Formula_1.2-4       UnitEvents_0.0.8   
[13] pec_2020.11.17      remotes_2.4.2       ggrepel_0.9.1       sessioninfo_1.1.1  
[17] numDeriv_2016.8-1.1 timereg_2.0.0       pillar_1.8.1        backports_1.2.1    
[21] lattice_0.20-44     glue_1.6.2          digest_0.6.29       gridtext_0.1.4     
[25] polyclip_1.10-0     checkmate_2.0.0     colorspace_2.0-3    Matrix_1.3-4       
[29] pkgconfig_2.0.3     devtools_2.4.2      purrr_0.3.4         mvtnorm_1.1-2      
[33] patchwork_1.1.1     scales_1.2.1        processx_3.5.2      tweenr_1.0.2       
[37] lava_1.6.9          ggforce_0.3.3       timechange_0.2.0    tibble_3.1.8       
[41] mgcv_1.8-36         generics_0.1.3      farver_2.1.0        usethis_2.1.6      
[45] ellipsis_0.3.2      cachem_1.0.5        withr_2.5.0         lazyeval_0.2.2     
[49] cli_3.6.0           survival_3.2-11     pammtools_0.5.7     magrittr_2.0.3     
[53] crayon_1.5.0        memoise_2.0.0       ggtext_0.1.1        ps_1.6.0           
[57] fs_1.5.0            fansi_1.0.4         nlme_3.1-152        MASS_7.3-54        
[61] xml2_1.3.2          pkgbuild_1.2.0      tools_4.1.3         prettyunits_1.1.1  
[65] lifecycle_1.0.3     stringr_1.5.0       munsell_0.5.0       callr_3.7.0        
[69] compiler_4.1.3      rlang_1.0.6         grid_4.1.3          iterators_1.0.13   
[73] rstudioapi_0.13     igraph_1.3.5        labeling_0.4.2      testthat_3.1.2     
[77] gtable_0.3.0        codetools_0.2-18    graphlayouts_0.7.1  R6_2.5.1           
[81] gridExtra_2.3       fastmap_1.1.0       utf8_1.2.3          rprojroot_2.0.2    
[85] desc_1.4.0          stringi_1.7.12      Rcpp_1.0.8          vctrs_0.5.2        
[89] tidyselect_1.2.0
```
