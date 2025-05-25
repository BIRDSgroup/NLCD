## How to process the raw yeast data  
Refer to `full_procedure.ipynb` in the `scripts` folder that uses the files present here to process the raw data. The processed yeast causal (`yeastgt_1_wilko1752_ready.txt`) and independent (`yeastgt_1_wilko1752_ready.txt`)data are present in the drive link. \
`mpmi_indices`: contains the indices of yeast trios (causal and independent) for 10 different seeds that were samples based on varying ranges of bcmi and spearman. created in `bcmi_regression_wilko.R`\
`var_indices` : contains the indices of yeast trios (causal and independent) for different unequal variances created in `yeast_variance.ipynb` \
`varandnonlinear_indices`: contains indices of trios that are intersecting the above two conditions (nonlinear and unvariance). Created in `yeast_variance.ipynb` \
`supp_causal_yeast_1752.csv` : contains the distance of all the causal yeast points from the regressed line between spearman and bcmi \
`supp_indep_yeast_1752.csv` : contains the distance of all the independent yeast points from the regressed line between spearman and bcmi 