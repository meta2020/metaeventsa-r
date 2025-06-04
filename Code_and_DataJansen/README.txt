TITLE OF THE MANUSCRIPT: Random-effects meta-analysis models for the odds ratio in the case of rare events under different data generating models: A simulation study

AUTHORS: Katrin Jansen and Heinz Holling

RESPONSIBLE FOR CODE: Katrin Jansen

CORRESPONDENCE TO: katrinjansen@uni-muenster.de
---------------------------------------------------------------------------------------------------------

1. FOLDER STRUCTURE:

The folder Code_and_Data has the following structure:

README.txt
/example
/manuscript
/simulation
    ./intermediate_results

The content of the subfolder simulation is further described in Section 2.
The content of the subfolder example is further described in Section 3.
The content of the subfolder manuscript is further described in Section 4.

2. SIMULATION:

The subfolder simulation contains the code required to reproduce the results of the simulation. This folder contains the scripts:

- 1_simordgm_simulation.R

- 2_simordgm_processing.R

- 3_simordgm_plots_manuscript.R

- 4_simordgm_plots_supp.R

- functions_simordgm.R

In section 2.2, we describe how to reproduce the simulation results by executing these scripts.

2.1 CONFIGURATIONS:

The scripts were tested under the following configurations:

- Simulation code tested under [scripts: 1_simordgm_simulation.R]:
 
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /Applic.HPC/Easybuild/skylake/2019a/software/MPI/intel/2019.1.144-GCC-8.2.0-2.31.1/impi/2018.4.274/R/3.6.0/lib64/R/lib/libR.so
LAPACK: /Applic.HPC/Easybuild/skylake/2019a/software/MPI/intel/2019.1.144-GCC-8.2.0-2.31.1/impi/2018.4.274/R/3.6.0/lib64/R/modules/lapack.so

locale:
[1] C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] numDeriv_2016.8-1.1 nleqslv_3.3.2       extraDistr_1.9.1   
 [4] BiasedUrn_1.07      metafor_3.0-2       lme4_1.1-21        
 [7] Matrix_1.2-17       dplyr_0.8.3         tidyr_1.0.0        
[10] doRNG_1.7.1         rngtools_1.3.1.1    pkgmaker_0.27      
[13] registry_0.5-1      doParallel_1.0.14   iterators_1.0.12   
[16] foreach_1.4.7      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3       mathjaxr_1.4-0   nloptr_1.2.1     pillar_1.4.2    
 [5] compiler_3.6.0   tools_3.6.0      boot_1.3-22      zeallot_0.1.0   
 [9] digest_0.6.23    nlme_3.1-140     tibble_2.1.3     lifecycle_0.1.0 
[13] lattice_0.20-38  pkgconfig_2.0.3  rlang_0.4.2      bibtex_0.4.2    
[17] withr_2.1.2      stringr_1.4.0    vctrs_0.2.0      grid_3.6.0      
[21] tidyselect_0.2.5 glue_1.3.1       R6_2.4.1         minqa_1.2.4     
[25] purrr_0.3.3      magrittr_1.5     backports_1.1.5  codetools_0.2-16
[29] splines_3.6.0    MASS_7.3-51.4    assertthat_0.2.1 xtable_1.8-4    
[33] stringi_1.4.3    crayon_1.3.4    
   
  
- Code for 2_simordgmprocessing.R tested under:

R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /Applic.HPC/Easybuild/skylake/2019a/software/MPI/intel/2019.1.144-GCC-8.2.0-2.31.1/impi/2018.4.274/R/3.6.0/lib64/R/lib/libR.so
LAPACK: /Applic.HPC/Easybuild/skylake/2019a/software/MPI/intel/2019.1.144-GCC-8.2.0-2.31.1/impi/2018.4.274/R/3.6.0/lib64/R/modules/lapack.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.0.0   dplyr_0.8.3   lme4_1.1-21   Matrix_1.2-17

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3       magrittr_1.5     splines_3.6.0    MASS_7.3-51.4   
 [5] tidyselect_0.2.5 lattice_0.20-38  R6_2.4.1         rlang_0.4.2     
 [9] minqa_1.2.4      grid_3.6.0       nlme_3.1-140     assertthat_0.2.1
[13] lifecycle_0.1.0  tibble_2.1.3     crayon_1.3.4     purrr_0.3.3     
[17] nloptr_1.2.1     vctrs_0.2.0      zeallot_0.1.0    glue_1.3.1      
[21] compiler_3.6.0   pillar_1.4.2     backports_1.1.5  boot_1.3-22     
[25] pkgconfig_2.0.3 


- Code for 3_simordgm_plots_manuscript.R tested under:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.2.0    ggthemes_4.2.4  ggh4x_0.2.1     forcats_0.5.1   stringr_1.4.0  
 [6] dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
[11] ggplot2_3.3.5   tidyverse_1.3.1

loaded via a namespace (and not attached):
 [1] cellranger_1.1.0 pillar_1.7.0     compiler_4.1.2   dbplyr_2.1.1     tools_4.1.2     
 [6] jsonlite_1.8.0   lubridate_1.8.0  lifecycle_1.0.1  gtable_0.3.0     pkgconfig_2.0.3 
[11] rlang_1.0.2      reprex_2.0.1     rstudioapi_0.13  DBI_1.1.2        cli_3.2.0       
[16] haven_2.5.0      xml2_1.3.3       withr_2.5.0      httr_1.4.2       fs_1.5.2        
[21] generics_0.1.2   vctrs_0.3.8      hms_1.1.1        grid_4.1.2       tidyselect_1.1.2
[26] glue_1.6.2       R6_2.5.1         fansi_1.0.2      readxl_1.4.0     tzdb_0.3.0      
[31] modelr_0.1.8     magrittr_2.0.1   backports_1.4.1  ellipsis_0.3.2   rvest_1.0.2     
[36] assertthat_0.2.1 colorspace_2.0-3 utf8_1.2.2       stringi_1.7.6    munsell_0.5.0   
[41] broom_0.8.0      crayon_1.5.1   

- Code for 4_simordgm_plots_manuscript.R tested under:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.2.0    ggthemes_4.2.4  ggh4x_0.2.1     forcats_0.5.1   stringr_1.4.0  
 [6] dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
[11] ggplot2_3.3.5   tidyverse_1.3.1

loaded via a namespace (and not attached):
 [1] cellranger_1.1.0 pillar_1.7.0     compiler_4.1.2   dbplyr_2.1.1     tools_4.1.2     
 [6] jsonlite_1.8.0   lubridate_1.8.0  lifecycle_1.0.1  gtable_0.3.0     pkgconfig_2.0.3 
[11] rlang_1.0.2      reprex_2.0.1     rstudioapi_0.13  DBI_1.1.2        cli_3.2.0       
[16] haven_2.5.0      xml2_1.3.3       withr_2.5.0      httr_1.4.2       fs_1.5.2        
[21] generics_0.1.2   vctrs_0.3.8      hms_1.1.1        grid_4.1.2       tidyselect_1.1.2
[26] glue_1.6.2       R6_2.5.1         fansi_1.0.2      readxl_1.4.0     tzdb_0.3.0      
[31] modelr_0.1.8     magrittr_2.0.1   backports_1.4.1  ellipsis_0.3.2   rvest_1.0.2     
[36] assertthat_0.2.1 colorspace_2.0-3 utf8_1.2.2       stringi_1.7.6    munsell_0.5.0   
[41] broom_0.8.0      crayon_1.5.1  

  
REQUIRED R PACKAGES:
 - Simulation [1_simordgm_simulation.R]:
	- doParallel
	- doRNG
	- tidyr
	- dplyr
	- lme4
	- metafor
	- BiasedUrn
	- extraDistr
	- nleqslv

 - Data processing [2_simordgm_processing.R]:
	- lme4
	- dplyr
	- tidyr

 - Plots for manucript [3_simordgm_plots_manuscript.R]:
	- tidyverse
	- ggh4x
	- ggthemes
	- scales

 - Plots for supplement [4_simordgm_plots_supp.R]:
	- tidyverse
	- ggh4x
	- ggthemes
	- scales


2.2 REPRODUCING THE SIMULATION RESULTS
 
Please make sure the required R packages are installed prior to running any of the scripts, since the scripts do not include code to 
install the packages on your computer. 

Below, we list each script in the order they should be run to conduct the simulation, process the simulation results, as well as to produce the figures.

The order in which the scripts should be run is also indicated by the numbers at the beginning of the respective file names. 
The folder also contains a file called "functions_simordgm.R", which contains functions needed in the script "1_simordgm_simulation.R". 
This file does not have to be run seperately but is sourced within the script "1_simordgm_simulation.R". 

The folder "simulation" also contains a subfolder called "intermediate_results" in which we have provided the results of the script "2_simordgm_processing.R". 
The results obtained from "1_simor_simulationdgm.R", which contain the raw simulation results, were too large to be uploaded. 

  - 1_simordgm_simulation.R : For this script, it is assumed that the working directory is set to the folder "simulation".
This file contains the code for the simulation. It returns 1152 rds files with the results for each of the 1152 simulation conditions, named "simor_X.rds" (X = 1, ..., 1152).
Per default, these are saved in a folder named "results" which will be created in the "simulation" folder if the script is run. 
The script prints statements about which condition is currently being computed. 
The script is written to be run on multiple cores in parallel and it will take several weeks to run (estimate of run time on computing cluster PALMA II of University of Muenster). 
Please also note that the .rds-files occupy a large amount of storage, which is also the reason why they are currently not provided on any public repository.
In case you want to run the script on your personal computer, please change the statement in line 33, cores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")), to
cores <- detectCores()-1. Then the script is run in parallel on all but one cores of your computer. In case you want to test the script with a smaller number of simulation replications,
please modify the statement in line 54, n_simulations <- 1000, and replace 1000 by the number of simulation replications you desire. We do not recommend to run the script on your
personal computer with more than three simulation replications.

  - 2_simordgm_processing.R : This file contains the code which prepares the simulation results so that they result in a smaller more manageable object with which we can then 
continue on to create the figures of the paper. The script may take several hours to run on your computer, depending on your device. We recommend not running the script on your personal computer. 
If you have run the script "1_simordgm_simulation.R", please set the working directory to the "results" folder which after having run the script "1_simordgm_simulation.R", 
contains the files "simor_X.rds" (for X = 1, ..., 1152). Furthermore, the file "design_simordgm.RData" is needed, which is produced by the script "1_simordgm_simulation.R" and 
consists of a data frame with the simulation conditions. Please make sure this file exists in the "results" folder. The script produces files named "results_X.RData" (X = 1, ..., 1152), 
which contain a data frame with the extracted results of the respective simulation condition. Finally, the script produces a summary of the "results_X.RData" files, 
which contains the performance measures for all simulation conditions and is named "simor_summary.RData".
The files named "results_X.RData" and "simor_summary.RData", which are produced by this script, can be downloaded from https://osf.io/xdthu/.
We also provide the files "simor_summary.RData" and "design_simordgm.RData" in the "intermediate_results" folder, since they are required to reproduce the figures.

  - 3_simordgm_plots_manuscript.R : This file contains the code that produces the figures shown in the manuscript. If you have run the scripts "1_simordgm_simulation.R" AND
"2_simordgm_processing.R" prior to executing this script, please set the working directory is set to the "results" folder. If you haven't run these scripts, please set
the working directory to the "intermediate_results" folder instead, which contains all required files.
Before figures are produced, the files "simor_summary.RData" and "design_simordgm.RData" are loaded, thus please make sure they exist in the working directory. 
Plots are saved to a folder names "plots_manuscript" which is created by the script in the working directory.
All plots are saved as eps files with their names indicating the figure number of the respective figure in the manuscript.

- 4_simordgm_plots_supp.R : This file contains the code that produces the plots shown in the supplement. If you have run the scripts "1_simordgm_simulation.R" AND
"2_simordgm_processing.R" prior to executing this script, please set the working directory is set to the "results" folder. If you haven't run these scripts, please set
the working directory to the "intermediate_results" folder instead, which contains all required files.
Before figures are produced, the files "simor_summary.RData" and "design_simordgm.RData" are loaded, thus please make sure they exist in the working directory. 
Plots are saved to a folder names "plots_supp" which is created by the script in the working directory (or already exists if you have run script "3_simordgm_plots_manuscript.R").
All plots are saved as pdf files with their names indicating the results they contain: 

Plots of descriptives are named 
- "avgnevents_kX.pdf" (plots of the average total number of events per arm), 
- "avg00_kX.pdf" (plots of the average number of single- and double-zero studies)
- "zeroarm_kX.pdf" (plots of the number of simulation runs in which either study arm or both study arms were zero)
where X is replaced by the number of studies which the results are depicted for. 

Plots of performance criteria are named
"criterion_grA_pcB_kC_nD_nE.pdf"
where "criterion" is replaced by the respective performance criterion (estimability, meanbias, medianbias, mse, coverage, ciwidth, mcse), A is replaced by the randomization ratio 
for which the results are depicted (1, 2), B is replaced by the event probability in the control group for which the results are depicted (01, 05),
C is replaced by the number of studies for which the results are depicted (5, 10, 30), D is replaced by the first sample size for which the results are depicted (50, 100),
and E is replaced by the second sample size for which the results are depicted (200, 500). 


3. EXAMPLE

The subfolder "example" contains the following files:

- 5_simor_example.R

- data_thomas.Rdata

- functions_example.R

- tablesB1andB2.R

These files are required to reproduce the analysis of the example (Section 6 of the manuscript).

The file "data_thomas.Rdata" contains data obtained from:

- Thomas, K. H., Martin, R. M., Knipe, D. W., Higgins, J. P. T., and Gunnel, D. (2015). Risk of neuropsychiatric
  adverse events associated with varenicline: systematic review and meta-analysis. BMJ 2015, 350:h1109. doi:
  10.1136/bmj.h1109

The data were obtained from the paper and the online supplement to the paper. These data are analysed in the script "5_simor_example.R". 
In section 3.2, it will be explained how the results for the example can be reproduced using this script.

The file "functions_example.R" does not have to be run. It is a file that contains functions required by the script "5_simor_example.R",
and is sourced in this script.

The script "tablesB1andB2.R" reproduces tables B1 and B2 from the manuscript. To execute this script, it is required that the 
working directory is set to the "example" subfolder, and that this subfolder contains the file "data_thomas.Rdata". 
This script requires that the xtable package is installed. Please make sure you have installed the package before running the script.

3.1 CONFIGURATIONS

- Code for 5_simor_example.R tested under:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] xtable_1.8-4        forcats_0.5.1       stringr_1.4.0       dplyr_1.0.8        
 [5] purrr_0.3.4         readr_2.1.2         tidyr_1.2.0         tibble_3.1.6       
 [9] ggplot2_3.3.5       tidyverse_1.3.1     lme4_1.1-28         nleqslv_3.3.2      
[13] numDeriv_2016.8-1.1 extraDistr_1.9.1    metafor_3.0-2       Matrix_1.3-4       
[17] devtools_2.4.3      usethis_2.1.5      

loaded via a namespace (and not attached):
 [1] httr_1.4.2        pkgload_1.2.4     jsonlite_1.8.0    splines_4.1.2     modelr_0.1.8     
 [6] brio_1.1.3        assertthat_0.2.1  cellranger_1.1.0  remotes_2.4.2     sessioninfo_1.2.2
[11] pillar_1.7.0      backports_1.4.1   lattice_0.20-45   glue_1.6.2        rvest_1.0.2      
[16] minqa_1.2.4       colorspace_2.0-3  pkgconfig_2.0.3   broom_0.8.0       haven_2.5.0      
[21] scales_1.2.0      processx_3.5.3    tzdb_0.3.0        generics_0.1.2    ellipsis_0.3.2   
[26] cachem_1.0.6      withr_2.5.0       cli_3.2.0         magrittr_2.0.1    crayon_1.5.1     
[31] readxl_1.4.0      memoise_2.0.1     ps_1.7.0          fs_1.5.2          fansi_1.0.2      
[36] nlme_3.1-153      MASS_7.3-54       xml2_1.3.3        pkgbuild_1.3.1    tools_4.1.2      
[41] prettyunits_1.1.1 hms_1.1.1         lifecycle_1.0.1   munsell_0.5.0     reprex_2.0.1     
[46] callr_3.7.0       compiler_4.1.2    rlang_1.0.2       grid_4.1.2        nloptr_2.0.0     
[51] rstudioapi_0.13   boot_1.3-28       testthat_3.1.3    gtable_0.3.0      DBI_1.1.2        
[56] R6_2.5.1          lubridate_1.8.0   fastmap_1.1.0     utf8_1.2.2        mathjaxr_1.6-0   
[61] rprojroot_2.0.3   desc_1.4.1        stringi_1.7.6     BiasedUrn_1.07    Rcpp_1.0.8       
[66] vctrs_0.3.8       dbplyr_2.1.1      tidyselect_1.1.2 


- Code for tablesB1andB2.R tested under:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] xtable_1.8-4

loaded via a namespace (and not attached):
[1] compiler_4.1.2 tools_4.1.2   


3.2 REPRODUCING THE EXAMPLE RESULTS

If not installed previously, you are required to install the following packages before you can run the script "5_simor_example.R":

- lme4
- metafor
- tidyverse
- xtable
- extraDistr (required by functions_example.R)
- nleqslv (required by functions_example.R)

In case you haven't installed these packages, please do so using the install.packages() command.

To reproduce the example results using the script "5_simor_example.R", make sure that the "example" subfolder contains the files "data_thomas.RData" and
"functions_example.R". Set the working directory to the "example" subfolder using the setwd() command in line 10 of the script.
and execute the script. The script extracts the data for the specific outcomes analysed in Section 6 (suicidal ideation and aggression), conducts all analyses and produces
the tables which contain the results. The (latex) tables are then saved in the files "tab_thomas2015_Suicidal Ideation.tex" (Table 4 of the manuscript)
and "tab_thomas2015_Aggression.tex" (Table 5 of the manuscript).
Throughout the process, the script produces a few warnings, which can be ignored.


4. MANUSCRIPT

The subfolder "manuscript" contains a single R script which reproduces Table 3 of the manuscript. The file which contains the R script is

- table3.R

This script requires that the xtable package is installed. Please make sure you have installed the package before running the script.

The script has been tested under the following configurations:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] xtable_1.8-4

loaded via a namespace (and not attached):
[1] compiler_4.1.2 tools_4.1.2   

