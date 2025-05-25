
The scripts folder contains :  \
-> `nlcd_main.py`: This file has all the functions defined to run NLCD \ 
-> `nlcd_simulation.py` : This file is used to run the simulation data and yeast data results  \
-> `nlcd_user.py`: This file imports the nlcd_main file and adds two more functions for the user, add seeds  \
-> `user_run.ipynb`: This notebook imports the above nlcd_user.py file and gives examples on how a user can run NLCD  \
-> `cit_run.r`: This file can be used to reproduce the results obtained from cit, can be used as a template to run the cit code \
-> `full_procedure.ipynb`: Code used to process the raw yeast data \
-> `yeast_variance.ipynb`: Contains code to perform liliferos test for variance, to generate indices based on variance and take the datasets based on variance and non linearity (bcmi) \
-> `extendofnonlinear.R` : Code to run spearman, correlation, mi and bcmi on yeast causal and independent dataset. Results present in the `yeast` folder: `mpmicorspear_causal_wilko1752.csv` and `mpmicorspear_causal_wilko1752.csv` \
-> `prepare_final_table.R`: Code used to find the causal trios, includes mappability mapping, p value adjustment and other filters applied \
-> `prepare_pvalues_epept_normal.R` : Compute the EPEPT final p-values from the epept intermediate values/losses \
-> `process_trio_automate.R`: Code used to compute human trissues from tissues using genotype and gene expression, includes using MatrixEQTL and covariate adjustment \
-> `bcmi_regression_wilko_git.R`: contains the code that plots the causal and independent data of yeast bcmi vs spearman CIT and NLCD. Also contains the code to select the indices for 10 seeds for a range of values values of bcmi and spearman, also includes the code for plotting Figure 4a and S8  \
-> `cit_code.R` : contains the code used to run CIT on simulation and yeast data \
Note: expression files are present in the drive link
