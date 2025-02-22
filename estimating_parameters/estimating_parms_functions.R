############################## Libraries ##############################
library(dplyr)

library(parallel)
devtools::load_all('/home/anjawu/thesis_work/hBayesDM/R') # on server
library(hBayesDM)
library(rstan)
library(tidyverse)
library(readxl)
library(future)
library(furrr)
library(ggpubr)
library(loo)


################################# Functions ####################################

################ MCMC Functions ################ 
# Function to extract numeric identifiers from the file name
## e.g. results_individuals_1-2-3-4-5-6-7.txt returns 1234567
extract_individuals <- function(file_name) {
  stringr::str_extract(file_name, "(?<=results_individuals_)(\\d+(?:-\\d+)*)") %>%
    stringr::str_replace_all("-", "")
}

#### Function to run the MCMC process
## Require: 
      # file_path: pathway to code and contained folders
      # data_files_path: pathway of where the data of individual IGT results are stored
      # file_name: name of text file that contains all individual IGT results to be run through MCMC
run_MCMC_function <- function(file_path, data_files_path, file_name) {
  print("Starting the MCMC function")
  datafile = paste0(data_files_path, file_name)
  print(paste("Processing file:", file_name))
  
  # Check if there is an MCMC_output folder to store the runs
  mcmc_output_path <- paste0(file_path, "MCMC_output/") 
  if (!dir.exists(mcmc_output_path)) {
    print("MCMC_output folder does not exist. Creating...")
    dir.create(mcmc_output_path)
    print("MCMC_output folder created.")
  }

  # remove the individuals which are being run in this MCMC run. 
  # This is specific for the file names in format: results_individuals_1-2-3-4-5-6-7.txt
  individuals <- extract_individuals(file_name)
  
  # number of iterations
  N <- 5000 
  # number of chains
  nCH <- 1
  for (runs in 1:5) {
    tryCatch({
      for (model in c("delta", "vpp", "orl")) {
        print(model)
        print(Sys.Date())
        T1 <- Sys.time()
        
        # Switch() https://www.codecademy.com/resources/docs/r/switch is an if/else statement that goes through three options.
        VBout <- switch(
          # model is input into the switch function, then depending on the three options, 
          # it runs the sigma inverse gamma MCMC hBayesDM (created by Anja Wu)
          model,
          "vpp" = igt_vpp_sigmainvgam(data = datafile, niter = N,
                                      nwarmup = N/4, nchain = nCH,
                                      ncore = nCH, nthin = 1, inits = "vb",
                                      indPars = "mean", modelRegressor = FALSE,
                                      vb = FALSE, adapt_delta = 0.95, stepsize = 1, 
                                      max_treedepth = 10),
          "orl" = igt_orl_sigmainvgam(data =  datafile, niter = N,
                                      nwarmup = N/4, nchain = nCH,
                                      ncore = nCH, nthin = 1, inits = "vb",
                                      indPars = "mean", modelRegressor = FALSE,
                                      vb = FALSE, adapt_delta = 0.95, stepsize = 1,
                                      max_treedepth = 10),
          igt_pvl_delta_sigmainvgam(data =  datafile, niter = N,
                                    nwarmup = N/4, nchain = nCH,
                                    ncore = nCH, nthin = 1, inits = "vb",
                                    indPars = "mean", modelRegressor = FALSE,
                                    vb = FALSE, adapt_delta = 0.95, stepsize = 1,
                                    max_treedepth = 10)
        )
        print("done. Elapsed time:")
        print(Sys.time() - T1)
        
        # Updated saving to include individuals (by number) within the name (as opposed to saying everyone like before)
        formatted_time <- format(Sys.time(), "%Y-%m-%d_%H-%M")
        saveRDS(VBout, file = paste0(mcmc_output_path, "igt_", individuals, "_", model, "_VBinit_MCMCrun_out_", formatted_time, ".rds"))
      }
    }, error = function(e) {
      print("Failed during processing.")
    })
  }
}

#### Function to pull best MCMC runs after MCMC is completed.
## This function goes through all combinations of individuals and their files
MCMC_best_pull_function <- function(file_path, num_of_indiv) {
  mcmc_output_path <- paste0(file_path, "MCMC_output/")
  all_files = list.files(mcmc_output_path)
  
  # Check if there is an MCMC_output/summary folder to store the runs
  mcmc_summary_path <- paste0(mcmc_output_path, "summary/")
  if (!dir.exists(mcmc_summary_path)) {
    print("MCMC_output/summary folder does not exist. Creating...")
    dir.create(mcmc_summary_path)
    print("MCMC_output/summary folder created.")
  }
  
  participant_numbers <- 1:num_of_indiv
  
  combinations <- unlist(lapply(1:num_of_indiv, function(n) combn(participant_numbers, n, simplify = FALSE)), recursive = FALSE)
  
  # Loop over each combination and to get best results for each combination of participants
  for (combo in combinations[8:length(combinations)]) { # must skip the solo numbers in the list of combos
    # this is where we look at each file with the specific number of individuals
    # parse through all_files for the ones that have the specific combo numbers
    individuals <- paste(combo, collapse = "")
    list_of_index_matches <- grep(paste0("igt_", individuals,"_"), all_files)
    # list of files corresponding to the specific combination of individuals
    f_use <- all_files[list_of_index_matches]
    
    gathered = NULL
    for(file_index in f_use){
      tryCatch({
        print(paste0("Using: ", file_index))
        file_in_use = readRDS(paste0(mcmc_output_path, file_index))
        
        model           = file_in_use$model
        
        #check if both chains are the same-ish  
        extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
        if(dim(extracted_samples)[2]==1){
          both_same_ish = NA
        } else{
          summaries_lok_post = apply(extracted_samples[,,"lp__"],2,summary)
          both_same_ish      = all.equal(summaries_lok_post["Median",1],summaries_lok_post["Median",2], tolerance = 100)
        }
        
        # combine and summarize
        lp_summary         = summary(as.vector(extracted_samples[,,"lp__"]))
        names(lp_summary)  = make.names(paste("lp__",names(lp_summary)))
        # gather LOOIC info from 
        Niter_per_chain  = dim(extracted_samples)[1]
        Nchain           = dim(extracted_samples)[2]
        fit_vals         = extract_ic(file_in_use)[[1]][[1]][,1]
        all_together     = tibble(filename         = file_index,
                                  model            = model,
                                  elpd_loo         = fit_vals["elpd_loo"],
                                  looic            = fit_vals["looic"],
                                  both_same_ish    = both_same_ish,
                                  Niter            = Niter_per_chain,
                                  Nchain           = Nchain,
                                  lp__.Min.        = lp_summary["lp__.Min."][[1]], 
                                  lp__.1st.Qu.     = lp_summary["lp__.1st.Qu."][[1]],
                                  lp__.Median      = lp_summary["lp__.Median"][[1]],
                                  lp__.Mean        = lp_summary["lp__.Mean"][[1]],
                                  lp__.3rd.Qu.     = lp_summary["lp__.3rd.Qu."][[1]],
                                  lp__.Max.        = lp_summary["lp__.Max."][[1]]
        )
        gathered = bind_rows(gathered,all_together)
        saveRDS(gathered, file = paste0(summary_path, "summary_interim.rds"))
      }, error = function(e){print(paste0("not this one ", file_index))}
      )
    }
    
    print("done_gathered... saving....")
    saveRDS(gathered, file = paste0(mcmc_summary_path, "summary_of_ALL_igt_subgroup_", individuals, "_VBinit_MCMC.rds"))
    
    #### keep track of the 'better' ones from all MCMC runs
    output = gathered  |> group_by(model) |> 
      arrange(desc(lp__.Median)) |> 
      mutate(maxlp = max(lp__.Median), same_as_best = near(maxlp, lp__.Median, tol = abs(.01*maxlp))) |> 
      ungroup()
    #  Find the best ones:
    print("grouping, now binding")
    better = output |> filter(same_as_best == TRUE)
    worse = output  |> filter(same_as_best == FALSE)
    
    saveRDS(better, file = paste0(mcmc_summary_path, "summary_of_BEST_igt_subgroup_", individuals, "_VBinit_MCMC.rds"))
    saveRDS(worse, file = paste0(mcmc_summary_path, "summary_of_NOTBEST_igt_subgroup_", individuals, "_VBinit_MCMC.rds"))
    print("saved, better and worse")
    
    #### MERGE THE RESULTS FROM THE BEST RUNS ####
    modelsets  = better|> select(model)|> unique()
    # SKIPPING PRE_ALLOCATION
    MCMC_RUNS = NULL
    for(model_grab in 1:nrow(modelsets)){
      f_use      <- better |> 
        filter( model           == modelsets[model_grab,"model"][[1]]) |> pull(filename)
      print("MCMC runs loop ")
      combined_samples = NULL
      for(file_index in f_use){ # f_use[1:min(15, length(f_use))]
        # here is where we discard burn in and concatenate all the MCMC runs.
        file_in_use = readRDS(paste0(mcmc_output_path, file_index ))
        
        if(is.null(combined_samples)){# the first one is saved
          combined_samples = file_in_use
        }else{# if it already exists, concat
          for(ind in 1:length(file_in_use$parVals)){
            if(!is.matrix(file_in_use$parVals[[ind]])){# it's a vector
              combined_samples$parVals[[ind]] = c(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]]) 
            }else{#2-dims
              combined_samples$parVals[[ind]] = rbind(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]])
            }
          }# end for over parVals
          combined_samples$fit = sflist2stanfit(list(combined_samples$fit, file_in_use$fit))
        }# end else
        
        # combine them all together in this format as well...
        model           = file_in_use$model
        #check if both chains are the same-ish  
        extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
        if(dim(extracted_samples)[2]==1){
          MCMC_RUNS        =   bind_rows(MCMC_RUNS,
                                         data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
        }else{   # eventually allow an arbitrary number of chains, but for now, just be lazy
          MCMC_RUNS        =   bind_rows(MCMC_RUNS,
                                         data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]),
                                         data.frame(extracted_samples[,2,],chain = 2, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
        }
        
      }
      combined_samples |> saveRDS(file = paste0(mcmc_summary_path, "merged_output_from_BEST_igt_subgroup_", individuals, "_",
                                                modelsets[model_grab,"model"][[1]],
                                                "_top15.rds"))
      extract_ic(combined_samples, ic="looic") |> saveRDS(file = paste0(mcmc_summary_path, "merged_output_from_BEST_igt_subgroup_", individuals, "_",
                                                                        modelsets[model_grab,"model"][[1]],
                                                                        "_top15_justlooic.rds"))
    }
    saveRDS(MCMC_RUNS, file = paste0(mcmc_summary_path, "samples_from_BEST__igt_subgroup_", individuals, "_VBinit_MCMC_top15.rds"))
    print("samples saved")
    
    ### Creating plots to assess fits of MCMC runs
    # Check if there is an MCMC_output/plots folder to store the runs
    mcmc_plots_path <- paste0(mcmc_output_path, "plots/") 
    if (!dir.exists(mcmc_plots_path)) {
      print("MCMC_output/plots folder does not exist. Creating...")
      dir.create(mcmc_plots_path)
      print("MCMC_output/plots folder created.")
    }
    
    for(model_grab in 1:nrow(modelsets)){
      print("plots over modelsets")
      dat2use = MCMC_RUNS |> filter(model == modelsets[model_grab,"model"][[1]]) 
      dat2use |> ggplot() +
        geom_point(aes(y=lp__, colour = filename, x = 1:nrow(dat2use)))
      ggsave(filename = paste0(mcmc_plots_path, "log_posterior_model_", individuals, "_",modelsets[model_grab,"model"][[1]],".png"))
    }
    
    #### Checking LOOIC results of BEST ####
    LOOIC_keeper = NULL
    loo_list = list()
    for(model_grab in 1:nrow(modelsets)){
      f_use      = MCMC_RUNS |>
        filter( model           == modelsets[model_grab,"model"][[1]]) |>
        pull(filename)|>unique()
      print("LOOIC runs loop ")
      for(file_index in f_use){
        print(file_index)
        file_in_use = readRDS(paste0(mcmc_output_path, file_index))
        loo = extract_ic(file_in_use, ic="looic")$LOOIC
        loo_list[[length(loo_list)+1]] = list(extract = loo, 
                                              file    = file_index,
                                              model   = modelsets[model_grab,"model"][[1]])
        
        delta_baseline_looic = data.frame(filename  = file_index, 
                                          model     = modelsets[model_grab,"model"][[1]], 
                                          LOOIC_est = loo$estimates["looic","Estimate"],
                                          LOOIC_se  = loo$estimates["looic","SE"])
        LOOIC_keeper = rbind(LOOIC_keeper, delta_baseline_looic)
      }
    }
    # gather LOOIC info from 
    LOOIC_keeper |> saveRDS( file = paste0(mcmc_summary_path, "LOOIC_from_BEST__igt_subgroup_", individuals, "_VBinit_MCMC_top15.rds"))
    loo_list     |> saveRDS( file = paste0(mcmc_summary_path, "LOOICLIST_from_BEST__igt_subgroup_", individuals, "_VBinit_MCMC_top15.rds"))
    
    # All possible comparisons between VPP and ORL models.
    LOOIC_keeper |> ggplot()+
      geom_boxplot(aes(y=LOOIC_est, x = model)) +
      ggsave(filename = paste0(mcmc_plots_path, "looic_boxplots_by_subgroup", individuals, "_.png"))
  }
}

# This function should be used if you have only run one file in MCMC run
MCMC_best_pull_single_file_function <- function(file_path) {
  mcmc_output_path <- paste0(file_path, "MCMC_output/") #output_location
  
  # Check if there is an MCMC_output/summary folder to store the runs
  mcmc_summary_path <- paste0(mcmc_output_path, "summary/") #summary_path
  if (!dir.exists(mcmc_summary_path)) {
    print("MCMC_output/summary folder does not exist. Creating...")
    dir.create(mcmc_summary_path)
    print("MCMC_output/summary folder created.")
  }
  
  # create list for files to use: 
  f_use = list.files(mcmc_output_path)
  
  gathered = NULL
  for(file_index in f_use){
    tryCatch({
      print(paste0("Using: ", file_index))
      file_in_use = readRDS(paste0(mcmc_output_path, file_index))

      model           = file_in_use$model
      
      #check if both chains are the same-ish  
      extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
      if(dim(extracted_samples)[2]==1){
        both_same_ish = NA
      } else{
        summaries_lok_post = apply(extracted_samples[,,"lp__"],2,summary)
        both_same_ish      = all.equal(summaries_lok_post["Median",1],summaries_lok_post["Median",2], tolerance = 100)
      }
      
      # combine and summarize
      lp_summary         = summary(as.vector(extracted_samples[,,"lp__"]))
      names(lp_summary)  = make.names(paste("lp__",names(lp_summary)))
      # gather LOOIC info from 
      Niter_per_chain  = dim(extracted_samples)[1]
      Nchain           = dim(extracted_samples)[2]
      fit_vals         = extract_ic(file_in_use)[[1]][[1]][,1]
      all_together     = tibble(filename         = file_index,
                                model            = model,
                                elpd_loo         = fit_vals["elpd_loo"],
                                looic            = fit_vals["looic"],
                                both_same_ish    = both_same_ish,
                                Niter            = Niter_per_chain,
                                Nchain           = Nchain,
                                lp__.Min.        = lp_summary["lp__.Min."][[1]], 
                                lp__.1st.Qu.     = lp_summary["lp__.1st.Qu."][[1]],
                                lp__.Median      = lp_summary["lp__.Median"][[1]],
                                lp__.Mean        = lp_summary["lp__.Mean"][[1]],
                                lp__.3rd.Qu.     = lp_summary["lp__.3rd.Qu."][[1]],
                                lp__.Max.        = lp_summary["lp__.Max."][[1]]
                          )
      gathered = bind_rows(gathered,all_together)
      saveRDS(gathered, file = paste0(summary_path, "summary_interim.rds"))
    }, error = function(e){print(paste0("not this one ", file_index))}
    )
  }
  
  print("done_gathered... saving....")
  saveRDS(gathered, file = paste0(mcmc_summary_path, "summary_of_ALL_igt_subgroup_VBinit_MCMC.rds"))
  
  #### keep track of the 'better' ones from all MCMC runs
  output = gathered  |> group_by(model) |> 
    arrange(desc(lp__.Median)) |> 
    mutate(maxlp = max(lp__.Median), same_as_best = near(maxlp, lp__.Median, tol = abs(.01*maxlp))) |> 
    ungroup()
  #  Find the best ones:
  print("grouping, now binding")
  better = output |> filter(same_as_best == TRUE)
  worse = output  |> filter(same_as_best == FALSE)
  
  saveRDS(better, file = paste0(mcmc_summary_path, "summary_of_BEST_igt_subgroup_VBinit_MCMC.rds"))
  saveRDS(worse, file = paste0(mcmc_summary_path, "summary_of_NOTBEST__igt_subgroup_VBinit_MCMC.rds"))
  print("saved, better and worse")
  
  #### MERGE THE RESULTS FROM THE BEST RUNS ####
  modelsets  = better|> select(model)|> unique()
  # SKIPPING PRE_ALLOCATION
  MCMC_RUNS = NULL
  for(model_grab in 1:nrow(modelsets)){
    f_use      <- better |> 
      filter( model           == modelsets[model_grab,"model"][[1]]) |> pull(filename)
    print("MCMC runs loop ")
    combined_samples = NULL
    for(file_index in f_use){ # f_use[1:min(15, length(f_use))]
      # here is where we discard burn in and concatenate all the MCMC runs.
      file_in_use = readRDS(paste0(mcmc_output_path, file_index ))

      if(is.null(combined_samples)){# the first one is saved
        combined_samples = file_in_use
      }else{# if it already exists, concat
        for(ind in 1:length(file_in_use$parVals)){
          if(!is.matrix(file_in_use$parVals[[ind]])){# it's a vector
            combined_samples$parVals[[ind]] = c(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]]) 
          }else{#2-dims
            combined_samples$parVals[[ind]] = rbind(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]])
          }
        }# end for over parVals
        combined_samples$fit = sflist2stanfit(list(combined_samples$fit, file_in_use$fit))
      }# end else
      
      # combine them all together in this format as well...
      model           = file_in_use$model
      #check if both chains are the same-ish  
      extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
      if(dim(extracted_samples)[2]==1){
        MCMC_RUNS        =   bind_rows(MCMC_RUNS,
                                       data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
      }else{   # eventually allow an arbitrary number of chains, but for now, just be lazy
        MCMC_RUNS        =   bind_rows(MCMC_RUNS,
                                       data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]),
                                       data.frame(extracted_samples[,2,],chain = 2, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
      }
      
    }
    combined_samples |> saveRDS(file = paste0(mcmc_summary_path, "merged_output_from_BEST_igt_subgroup_",
                                              modelsets[model_grab,"model"][[1]],
                                              "_top15.rds"))
    extract_ic(combined_samples, ic="looic") |> saveRDS(file = paste0(mcmc_summary_path, "merged_output_from_BEST_igt_subgroup_",
                                                                      modelsets[model_grab,"model"][[1]],
                                                                      "_top15_justlooic.rds"))
  }
  saveRDS(MCMC_RUNS, file = paste0(mcmc_summary_path, "samples_from_BEST__igt_subgroup_VBinit_MCMC_top15.rds"))
  print("samples saved")
  
  ### Creating plots to assess fits of MCMC runs
  # Check if there is an MCMC_output/plots folder to store the runs
  mcmc_plots_path <- paste0(mcmc_output_path, "plots/") 
  if (!dir.exists(mcmc_plots_path)) {
    print("MCMC_output/plots folder does not exist. Creating...")
    dir.create(mcmc_plots_path)
    print("MCMC_output/plots folder created.")
  }
  
  for(model_grab in 1:nrow(modelsets)){
    print("plots over modelsets")
    dat2use = MCMC_RUNS |> filter(model == modelsets[model_grab,"model"][[1]]) 
    dat2use |> ggplot() +
                geom_point(aes(y=lp__, colour = filename, x = 1:nrow(dat2use)))
    ggsave(filename = paste0(mcmc_plots_path, "log_posterior_model_",modelsets[model_grab,"model"][[1]],".png"))
  }
  
  #### Checking LOOIC results of BEST ####
  LOOIC_keeper = NULL
  loo_list = list()
  for(model_grab in 1:nrow(modelsets)){
    f_use      = MCMC_RUNS |>
                  filter( model           == modelsets[model_grab,"model"][[1]]) |>
                    pull(filename)|>unique()
    print("LOOIC runs loop ")
    for(file_index in f_use){
      print(file_index)
      file_in_use = readRDS(paste0(mcmc_output_path, file_index))
      loo = extract_ic(file_in_use, ic="looic")$LOOIC
      loo_list[[length(loo_list)+1]] = list(extract = loo, 
                                            file    = file_index,
                                            model   = modelsets[model_grab,"model"][[1]])
      
      delta_baseline_looic = data.frame(filename  = file_index, 
                                        model     = modelsets[model_grab,"model"][[1]], 
                                        LOOIC_est = loo$estimates["looic","Estimate"],
                                        LOOIC_se  = loo$estimates["looic","SE"])
      LOOIC_keeper = rbind(LOOIC_keeper, delta_baseline_looic)
    }
  }
  # gather LOOIC info from 
  LOOIC_keeper |> saveRDS( file = paste0(mcmc_summary_path, "LOOIC_from_BEST__igt_subgroup_VBinit_MCMC_top15.rds"))
  loo_list     |> saveRDS( file = paste0(mcmc_summary_path, "LOOICLIST_from_BEST__igt_subgroup_VBinit_MCMC_top15.rds"))
  
  # All possible comparisons between VPP and ORL models.
  LOOIC_keeper |> ggplot()+
    geom_boxplot(aes(y=LOOIC_est, x = model)) +
    ggsave(filename = paste0(mcmc_plots_path, "looic_boxplots_by_subgroup.png"))
} 



################ Calculating Functions ################ 
# Function to calculate mean and variance for uniform distribution 
calc_uniform_mean <- function(a, b) {
  mean <- (a + b)/2
  return(mean)
}

calc_uniform_sd <- function(a, b) {
  var <- (b - a)^2 / 12
  sd <- sqrt(var)
  return(sd)
}

calculate_inverse_gamma_parameters <- function(sigma) {
  a <- (2*pi -2)/(pi -2)
  b <- sqrt((a-1)^2*(a-2)*sigma^2*(1-2/pi))
  parameters <- list('a' = a, 'b' = b)
  return(parameters)
}

calc_mean_for_Sigma <- function(sigma) {
  parameters <- calculate_inverse_gamma_parameters(sigma)
  a <- parameters$a
  b <- parameters$b
  
  mean <- b/(a-1) # this returns the variance
  return(mean)
}

bic_indiv_calc_fn <- function(log_lik, n, n_parms) {
  individual = -2*log_lik
  penalty = n_parms*log(n)
  bic_individual = individual + penalty
  return(bic_individual)
}



################ Estimating Parameters Functions ################ 
## Collect all priors from all three models
model_priors_function <- function(model_results, half_n_sigma, half_c_sigma){
  # Initialize the list to store all data frames
  all_priors <- list()
  
  # VPP DataFrame https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_vpp.stan
  vpp <- data.frame(
    A = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b), 
    alpha = c(calc_uniform_mean(0, 2), calc_uniform_sd(0, 2),  calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    cons = c(calc_uniform_mean(0, 5), calc_uniform_sd(0, 5), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    lambda = c(calc_uniform_mean(0, 10), calc_uniform_sd(0, 10), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    epP = c(0, 1, calc_mean_for_Sigma(half_c_sigma), calculate_inverse_gamma_parameters(half_c_sigma)$b),
    epN = c(0, 1, calc_mean_for_Sigma(half_c_sigma), calculate_inverse_gamma_parameters(half_c_sigma)$b),
    K = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    w = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b)
  )
  rownames(vpp) <- c("mean", "sd", "Sigma", "Sigma_b")
  model_results$vpp$prior <- vpp
  # Number of individuals is important, so store as variable
  model_results$vpp$total_individuals <- length(model_mcmc_results$vpp$allIndPars$subjID)
  
  # PVL_Delta DataFrame https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_pvl_delta.stan
  pvl_delta <- data.frame(
    A = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b), 
    alpha = c(calc_uniform_mean(0, 2), calc_uniform_sd(0, 2),  calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    cons = c(calc_uniform_mean(0, 5), calc_uniform_sd(0, 5), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    lambda = c(calc_uniform_mean(0, 10), calc_uniform_sd(0, 10), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b)
  )
  rownames(pvl_delta) <- c("mean", "sd", "Sigma", "Sigma_b")
  model_results$pvl_delta$prior <- pvl_delta
  model_results$pvl_delta$total_individuals <- length(model_mcmc_results$pvl_delta$allIndPars$subjID)
  
  # ORL DataFrame https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_orl.stan
  orl <- data.frame(
    Arew = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    Apun = c(calc_uniform_mean(0, 1), calc_uniform_sd(0, 1), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    K = c(calc_uniform_mean(0, 5), calc_uniform_sd(0, 5), calc_mean_for_Sigma(half_n_sigma), calculate_inverse_gamma_parameters(half_n_sigma)$b),
    betaF = c(0, 1, calc_mean_for_Sigma(half_c_sigma), calculate_inverse_gamma_parameters(half_c_sigma)$b),
    betaP = c(0, 1, calc_mean_for_Sigma(half_c_sigma), calculate_inverse_gamma_parameters(half_c_sigma)$b)
  )
  rownames(orl) <- c("mean", "sd", "Sigma", "Sigma_b")
  model_results$orl$prior<- orl
  model_results$orl$total_individuals <- length(model_mcmc_results$orl$allIndPars$subjID)
  
  return(model_results)
}


## Collecting all information from the mcmc results and adding to the model_results df
mcmc_results_amalgamation_function <- function(model_mcmc_results, model_results) {
  # Cycle through each model
  for (model in names(model_results)) {
    # Update name to be including _top15, as this is how it is saved in our data
    model_name <- paste0(model, "_top15")
    
    # Get a parameters list based on the model examined; exclude the subjID
    parameter_names_list <- colnames(model_mcmc_results[[model_name]]$allIndPars)[2:ncol(model_mcmc_results[[model_name]]$allIndPars)]
    # this summary accesses the mean and standard deviation for each parameter both population and individual parameters
    fit_summary <- rstan::summary(model_mcmc_results[[model_name]]$fit)$summary
    
    # Create a data frame to store population and individual statistics
    population_results_df <- data.frame(matrix(ncol = length(parameter_names_list), nrow = 2))
    individual_mean_df <- data.frame(matrix(ncol = length(parameter_names_list), nrow = 0))
    individual_sd_df <- data.frame(matrix(ncol = length(parameter_names_list), nrow = 0))
    loglik_df <- data.frame(matrix(ncol = 1, nrow = 0))
    
    # populate the population mean for each parameter from fit summary
    population_mean_row <- data.frame(t(sapply(parameter_names_list, function(parameter) {
      fit_summary[paste0('mu_', parameter), 'mean']
    })))
    
    # populate the population standard deviation for each parameter from fit summary
    population_std_row <- data.frame(t(sapply(parameter_names_list, function(parameter) {
      fit_summary[paste0('mu_', parameter), 'sd']
    })))
    # must update column names as they aren't saved as parameters
    colnames(population_std_row) <- parameter_names_list
    
    # populate the population standard deviation for each parameter from fit summary
    population_Sigma_row <- data.frame(t(sapply(seq_along(parameter_names_list), function(p) {
      fit_summary[paste0('sigma[', p, ']'), 'mean']
    })))
    # must update column names as they aren't saved as parameters
    colnames(population_Sigma_row) <- parameter_names_list
    
    # combining the mean and standard deviation, we make one dataframe for our population parameters:
    population_results_df <- rbind(population_mean_row, population_std_row, population_Sigma_row)
    rownames(population_results_df)[1:3] <- c("mean", "sd", "Sigma")
    
    # Now we cycle through each individual
    for (i in 1:model_results[[model]]$total_individuals) {
      # collect each individual mean
      individual_mean_row <- data.frame(t(sapply(parameter_names_list, function(parameter) {
        # Select individual and parameter, as this is how it is stored in the fit summary
        individual_param_selected <- paste0(parameter,'[',i,']')
        fit_summary[individual_param_selected,'mean']
      })))
      # combine the individual mean row by row
      individual_mean_df <- rbind(individual_mean_df, individual_mean_row)
      rownames(individual_mean_df)[i] <- i
      
      # collect each individual standard deviation
      individual_sd_row <- data.frame(t(sapply(parameter_names_list, function(parameter) {
        # Select individual and parameter, as this is how it is stored in the fit summary
        individual_param_selected <- paste0(parameter,'[',i,']')
        fit_summary[individual_param_selected,'sd']
      })))
      # combine the individual standard deviation row by row
      individual_sd_df <- rbind(individual_sd_df, individual_sd_row)
      rownames(individual_sd_df)[i] <- i
      
      # Create new dataframe to store log-likelihood and bic score
      individual_log_lik <- fit_summary[paste0('log_lik[',i,']'),'mean']
      # use function to calculate each bic score
      individual_bic <- bic_indiv_calc_fn(individual_log_lik, 100, length(parameter_names_list))
      log_bic_column <- cbind(individual_log_lik, individual_bic)
      colnames(log_bic_column) <- c("log_lik", "bic")
      loglik_df <- rbind(loglik_df, log_bic_column)
      rownames(loglik_df)[i] <- i
    }
    
    # store each model results into it's own object:
    model_results[[model]]$population_results <-  population_results_df
    model_results[[model]]$individual$mean <-  individual_mean_df
    model_results[[model]]$individual$sd <-  individual_sd_df
    model_results[[model]]$individual$criterion <-  loglik_df
  }
  return(model_results)
}


## This returns the initial results, but with the updated population parameters
population_updating_fn <- function(starting_model_results, best_fit_set) {
  # copy original model results so we can return updated population version
  population_update <- starting_model_results
  
  # cycle through model
  for (model in names(starting_model_results)) {
    model_selected <- starting_model_results[[model]]
    m <- length(best_fit_set[[model]])
    # update number of individuals in each model
    population_update[[model]]$total_individuals <- m
    
    # Create temporary dataframe that stores each models updated parameters
    model_updates <- matrix(nrow = 3, ncol = model_selected$parameter_count)
    colnames(model_updates) <- colnames(model_selected$prior)
    rownames(model_updates) <- c('mean', 'sd', 'Sigma')
    model_updates <- data.frame(model_updates)
    
    # cycle through parameters within model
    for (parameter in colnames(model_selected$prior)) {
      # need mu prior for each parameter mean distribution
      mu_omega0 <- model_selected$prior['mean', parameter]
      # need sigma prior for each parameter mean distribution
      sigma_omega0 <- model_selected$prior['sd', parameter]
      # need population variance initial value
      Sigma_omega <- model_selected$prior['Sigma', parameter]
      # Sigma_omega <- model_selected$population_results['Sigma', parameter]
      # mcmc run population mean:
      theta_omega <- model_selected$population_results['mean', parameter]
      
      ### Equation (19):
      ### a and b need to be updated to real values!!
      Sigma_a <- (2*pi -2)/(pi -2) # consistent
      Sigma_b <- model_selected$prior['Sigma_b', parameter] #parameter dependent
      
      # initiate sum for all individuals by parameter for whichever are within the specific model
      omega_population_sum <- 0
      # initiate sum for (theta_omega - omega_i)^2 for all i in m-set:
      theta_omega_population_summed <- 0
      # cycle through individuals within the model (from best_fit_set[[model]])
      for (individual in best_fit_set[[model]]) {
        # need all individual parameters added together
        omega_i <- model_selected$individual$mean[as.character(individual), parameter]
        omega_population_sum <- omega_population_sum + omega_i

        # calculating sum of (\theta_\omega - \omega_i)^2
        theta_omega_population <- (theta_omega - omega_i)^2
        theta_omega_population_summed <- theta_omega_population_summed + theta_omega_population
      }
      
      ### Equation (3):
      mean_numerator <- sigma_omega0^2*omega_population_sum + Sigma_omega*mu_omega0
      mean_denominator <- m*sigma_omega0^2 + Sigma_omega
      # update mean for the specific parameter in the temporary model dataframe
      model_updates['mean',parameter] <- mean_numerator/mean_denominator
      
      var_numerator <- Sigma_omega*sigma_omega0^2
      var_denominator <- m*sigma_omega0^2 + Sigma_omega
      # update variance for the specific parameter in the temporary model dataframe
      # NOTE: This is sqrt() as is written in equation (3)
      model_updates['sd',parameter] <- sqrt(var_numerator/var_denominator)
      
      # make equation (4)
      sigma_numerator <- theta_omega_population_summed/2 + Sigma_b
      sigma_denomiator <- m/2 + Sigma_a -1
      model_updates['Sigma',parameter] <- sigma_numerator/sigma_denomiator
    }
    population_update[[model]]$population_results <- model_updates
  }
  return(population_update)
}



## This returns the initial results, but with the updated population parameters and updated individual parameters
individual_updating_fn <- function(starting_model_results, best_fit_set, population_update) {
  # I want to copy the update population model results and then update the individual parameters as well
  individual_update <- population_update

  # cycle through model
  for (model in names(starting_model_results)) {
    model_selected <- starting_model_results[[model]]
    m <- length(best_fit_set[[model]])
    
    # Create temporary dataframe that stores each models updated parameters mean
    mean_model_updates <- matrix(nrow = m, ncol = model_selected$parameter_count)
    colnames(mean_model_updates) <- colnames(model_selected$prior)
    mean_model_updates <- data.frame(mean_model_updates)
    # Create temporary dataframe that stores each models updated parameters var
    sd_model_updates <- matrix(nrow = m, ncol = model_selected$parameter_count)
    colnames(sd_model_updates) <- colnames(model_selected$prior)
    sd_model_updates <- data.frame(sd_model_updates)
    
    # cycle through parameters within model
    for (parameter in colnames(model_selected$prior)) {
      # need original population mean value
      theta_omega <- model_selected$population_results['mean', parameter]
      # need updated population mean value:
      hat_theta_omega <- population_update[[model]]$population_results['mean', parameter]
      # population-level variance:
      Sigma_omega <- model_selected$prior['Sigma', parameter]
      # updated population variance
      hat_Sigma_omega <- population_update[[model]]$population_results['Sigma', parameter]
      
      # to keep adding row by row for mean and var in dataframe
      i = 1
      # cycle through individuals within the model (from best_fit_set[[model]])
      for (individual in best_fit_set[[model]]) {
        # need all individual parameters added together omega_i = mu_i??
        omega_i <- model_selected$individual$mean[as.character(individual), parameter]
        sigma_i <- model_selected$individual$sd[as.character(individual), parameter]
        
        # make equation (13)
        mean_numerator <- omega_i*hat_Sigma_omega*Sigma_omega + hat_theta_omega*Sigma_omega*sigma_i^2 - theta_omega*hat_Sigma_omega*sigma_i^2
        mean_denominator <- hat_Sigma_omega*Sigma_omega + sigma_i^2*Sigma_omega - sigma_i^2*hat_Sigma_omega
        # update mean for the specific parameter in the temporary model dataframe
        mean_model_updates[i, parameter] <- mean_numerator/mean_denominator
        
        var_numerator <- hat_Sigma_omega*Sigma_omega*sigma_i^2
        var_denominator <- hat_Sigma_omega*Sigma_omega + sigma_i^2*Sigma_omega - sigma_i^2*hat_Sigma_omega
        # update variance for the specific parameter in the temporary model dataframe
        # NOTE: This is sqrt() as is written in equation (13)
        sd_model_updates[i, parameter] <- sqrt(var_numerator/var_denominator)
        
        i = i+1
      }
      rownames(mean_model_updates) <- best_fit_set[[model]]
      rownames(sd_model_updates) <- best_fit_set[[model]]
    }
    individual_update[[model]]$individual$mean <- mean_model_updates
    individual_update[[model]]$individual$sd <- sd_model_updates # changed to sd to match what it was before
  }
  return(individual_update)
}







