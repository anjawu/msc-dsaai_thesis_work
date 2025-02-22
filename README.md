
# File breakdown with description

All the code in this repository accompanies the thesis paper: 

## Individual IGT Data Analysis
Contained within folder "individuals"

- combining_IGT_samples.Rmd: contains code to take a folder with .csv individual IGT results and create .txt files with each combination of individuals as their own separate .txt files. This is to be used for both approaches.


## Estimating Parameters Approach
Contained within folder "estimating_parameters"

-  estimating_parms_functions.R: contains all functions needed to run the algorithm #3 from the paper; updating placement based on estimating parameters
-  estimating_parms_run.Rmd: implements all functions to run MCMC, collect results and create the stabilizing loop that is at the crux of algorithm #3
-  individual_igt_w_inverseGamma/collect_individual_bics.Rmd: contains code that pulls all MCMC results for all combinations of individuals in the sample and creates box plots to show distributions of BIC score for each individual for the new added models to the hBayesDM (inverse gamma models)
	-  Results of box plots are stored within folder "individual_bic_scores/boxplots/"

**Please note** that the hBayesDM library that needs to be used for this section is my forked library: https://github.com/anjawu/hBayesDM, **not** the original hBayesDM.

## Variation Bayes Approach
Contained within folder "vb"

-  vb_functions.R: contains all functions needed to run the algorithm #4 from the paper; variation bayes method
-  vb_run.Rmd: implements all functions to run MCMC, collect results and create the stabilizing loop that is at the crux of algorithm #4
-  individual_igt/collect_individual_bics_vb.Rmd: contains code that pulls all MCMC results for all combinations of individuals in the sample and creates box plots to show distributions of BIC score for each individual for the original hBayesDM distribution.
	-  Results of box plots are stored within folder "individual_bic_scores/boxplots/"


