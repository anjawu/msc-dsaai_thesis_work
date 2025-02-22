############################## Libraries ##############################
library(dplyr)
library(parallel)
library(hBayesDM) # use their original library
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
          "vpp" = igt_vpp(data = datafile, niter = N,
                                      nwarmup = N/4, nchain = nCH,
                                      ncore = nCH, nthin = 1, inits = "vb",
                                      indPars = "mean", modelRegressor = FALSE,
                                      vb = FALSE, adapt_delta = 0.95, stepsize = 1, 
                                      max_treedepth = 10),
          "orl" = igt_orl(data =  datafile, niter = N,
                                      nwarmup = N/4, nchain = nCH,
                                      ncore = nCH, nthin = 1, inits = "vb",
                                      indPars = "mean", modelRegressor = FALSE,
                                      vb = FALSE, adapt_delta = 0.95, stepsize = 1,
                                      max_treedepth = 10),
          igt_pvl_delta(data =  datafile, niter = N,
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
# calculating individual BIC score
bic_indiv_calc_fn <- function(log_lik, n, n_parms) {
  individual = -2*log_lik
  penalty = n_parms*log(n)
  bic_individual = individual + penalty
  return(bic_individual)
}


# get best placement of BIC
get_best_placement_function <- function(mcmc_summary_path, individuals_string_with_space) {
  individuals_list = unlist(strsplit(individuals_string_with_space, split = " "))
  individuals_string = paste0(individuals_list, collapse = "")
  ### Read in all model MCMC outputs ###
  model_mcmc_results <- list()
  model_names <- c('vpp', 'orl', 'pvl_delta')#names(model_results)
  
  for (m in 1:length(model_names)) {
    file_to_import <- paste0(mcmc_summary_path, "merged_output_from_BEST_igt_subgroup_", individuals_string, "_igt_", model_names[m], "_top15.rds")
    print(paste0("Importing file: ", file_to_import))
    model_mcmc_results[[paste0(model_names[m], "_top15")]]  <- readRDS(file = file_to_import)
  }
  
  bic_results <- list()
  # Collect loglikelihoods and calculate BIC
  for (model in model_names) {
    # Update name to be including _top15, as this is how it is saved in our data
    model_name <- paste0(model, "_top15")
    
    # Get a parameters list based on the model examined; exclude the subjID
    parameter_names_list <- colnames(model_mcmc_results[[model_name]]$allIndPars)[2:ncol(model_mcmc_results[[model_name]]$allIndPars)]
    # this summary accesses the mean and standard deviation for each parameter both population and individual parameters
    fit_summary <- rstan::summary(model_mcmc_results[[model_name]]$fit)$summary
    
    temp_bic_df <- data.frame(matrix(ncol = 1, nrow = 0))
    # Now we cycle through each individual
    for (i in 1:nrow(model_mcmc_results[[model_name]]$allIndPars)) {
      # Create new dataframe to store log-likelihood and bic score
      individual_log_lik <- fit_summary[paste0('log_lik[',i,']'),'mean']
      # use function to calculate each bic score
      individual_bic <- bic_indiv_calc_fn(individual_log_lik, 100, length(parameter_names_list))
      # log_bic_column <- cbind(individual_log_lik, individual_bic)
      # colnames(log_bic_column) <- c("log_lik", "bic")
      # colnames(individual_bic) <- c("bic")
      temp_bic_df <- rbind(temp_bic_df, individual_bic)
      rownames(temp_bic_df)[i] <- i
    }
    colnames(temp_bic_df) <- "bic"
    bic_results[[model]] <- temp_bic_df
  }

  initial_best_fit_set <- list(vpp = c(), orl = c(), pvl_delta = c())
  
  for (indiv in 1:length(individuals_list)) {
    # compare columns: vpp_BIC, orl_BIC, pvl_delta_BIC for lowest value:
    vpp_bic <- bic_results$vpp['bic'][indiv,]
    orl_bic <- bic_results$orl['bic'][indiv,]
    pvl_delta_bic <- bic_results$pvl_delta['bic'][indiv,]
    
    ## find model with smallest bic:
    bic_min_index <- which.min(c(vpp_bic, orl_bic, pvl_delta_bic))
    ## find model name that has the lowest bic:
    lowest_bic_model_name <- model_names[bic_min_index]
    
    ## initial_best_fit_set$vpp # has list of all individual parameters that are the best for vpp
    ## initial_best_fit_set$orl # has list of all individual parameters that are the best for orl
    ## initial_best_fit_set$pvl_delta # has list of all individual parameters that are the best for pvl_delta
    initial_best_fit_set[[lowest_bic_model_name]] <- c(initial_best_fit_set[[lowest_bic_model_name]], individuals_list[indiv])
  }
  return(initial_best_fit_set)
}

################ VB Estimating Functions ################ 
## VB estimation runs
# This function is to use the variational inference to estimate a Normal from the data sample created.
# input only file name of results to run loop on
run_vb_estimate_function <- function(data_files_path, vb_file_path, file_name) {
  print(paste("Processing file:", file_name))
  
  datafile = paste0(data_files_path, file_name)
  print(paste("File path:", datafile))
  individuals <- extract_individuals(file_name)
  print(individuals)
  
  N <- 100
  nCH <- 1
  for (runs in 1:5) {
    tryCatch({
      print("STARTING")
      for (model in c("delta", "vpp", "orl")) {
        print(model)
        print(Sys.Date())
        T1 <- Sys.time()
        
        #  Here I was recommended the switch() https://www.codecademy.com/resources/docs/r/switch 
        VBout <- switch(
          model,
          "vpp" = igt_vpp(data = datafile, niter = N,
                          nwarmup = N/4, nchain = nCH,
                          ncore = nCH, nthin = 1, inits = "vb",
                          indPars = "mean", modelRegressor = FALSE,
                          vb = TRUE, adapt_delta = 0.95, stepsize = 1, #vb = TRUE
                          max_treedepth = 10),
          "orl" = igt_orl(data =  datafile, niter = N,
                          nwarmup = N/4, nchain = nCH,
                          ncore = nCH, nthin = 1, inits = "vb",
                          indPars = "mean", modelRegressor = FALSE,
                          vb = TRUE, adapt_delta = 0.95, stepsize = 1,
                          max_treedepth = 10),
          igt_pvl_delta(data =  datafile, niter = N,
                        nwarmup = N/4, nchain = nCH,
                        ncore = nCH, nthin = 1, inits = "vb",
                        indPars = "mean", modelRegressor = FALSE,
                        vb = TRUE, adapt_delta = 0.95, stepsize = 1,
                        max_treedepth = 10)
        )
        print("done. Elapsed time:")
        print(Sys.time() - T1)
        
        formatted_time <- format(Sys.time(), "%Y-%m-%d_%H-%M")
        # Updated saving to include individuals (by number) within the name (as opposed to saying everyone like before)
        saveRDS(VBout, file = paste0(vb_file_path, "vb_output/igt_", individuals, "_", model, "_VBinit_estimates_", formatted_time, ".rds"))
      }
    }, error = function(e) {
      print("Failed during processing.")
    })
  }
}

## VB estimation merging best results
# input individual combination that you want to specifically search VB results for 
### f_use <- only select specific with the individuals; figure out how to do this
# individuals = '1234567'
run_vb_merge_results_function <- function(vb_file_path, individuals) {
  # get correct VB outputs for individuals you are looking for:
  vb_output_path = paste0(vb_file_path, "vb_output/")
  vb_summary_path = paste0(vb_output_path, "summary/")
  all_files = list.files(vb_output_path)
  list_of_index_matches <- grep(paste0("igt_", individuals,"_"), all_files)
  # list of files corresponding to the specific combination of individuals
  f_use <- all_files[list_of_index_matches]
  
  gathered = NULL
  for(file_index in f_use){tryCatch({
    file_in_use = readRDS(paste0(vb_output_path, file_index))
    
    
    model           = file_in_use$model
    print(paste0("File: ", file_index, ", model: ", model))
    #check if both chains are the same-ish  
    extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
    # since this version of variational bayes does not return lp_approx, we need to estimate it by taking the sum of the log_lik
    # first collect all log-likelihoods from sample:
    log_lik_values <- rstan::extract(file_in_use$fit, pars = "log_lik")$log_lik
    # sum over each row (or sample)
    lp_approx_values <- rowSums(log_lik_values)
    # take summary to get final lp_approx_
    lp_summary <- summary(lp_approx_values)
    
    if(dim(extracted_samples)[2]==1){
      both_same_ish = NA
    }else{
      # not sure how to fix this, but it doesn't seem to reach it, so it is okay.
      summaries_lok_post = lp_summary#apply(extracted_samples[,,"lp__"],2,summary)
      both_same_ish      = all.equal(summaries_lok_post["Median",1],summaries_lok_post["Median",2], tolerance = 100)
    }
    
    # combine and summarize
    names(lp_summary)  = make.names(paste("lp_approx__",names(lp_summary)))
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
                              lp_approx__.Min.        = lp_summary["lp_approx__.Min."][[1]], 
                              lp_approx__.1st.Qu.     = lp_summary["lp_approx__.1st.Qu."][[1]],
                              lp_approx__.Median      = lp_summary["lp_approx__.Median"][[1]],
                              lp_approx__.Mean        = lp_summary["lp_approx__.Mean"][[1]],
                              lp_approx__.3rd.Qu.     = lp_summary["lp_approx__.3rd.Qu."][[1]],
                              lp_approx__.Max.        = lp_summary["lp_approx__.Max."][[1]]
    )
    gathered = bind_rows(gathered, all_together)
    saveRDS(gathered, file = paste0(vb_summary_path, "summary_", individuals, "_interim.rds"))
  }, error = function(e){print(paste0("not this one ", file_index))}) 
  }
  
  print("done_gathered... saving....")
  saveRDS(gathered, file = paste0(vb_summary_path, "summary_of_ALL_igt_subgroup_", individuals, "_VBinit_vb_estimate.rds"))
  
  gathered = readRDS(paste0(vb_summary_path, "summary_of_ALL_igt_subgroup_", individuals, "_VBinit_vb_estimate.rds"))
  
  #### keep track of the 'better' ones
  output = gathered  |> group_by(model) |> 
    arrange(desc(lp_approx__.Median)) |> 
    mutate(maxlp = max(lp_approx__.Median), same_as_best = near(maxlp, lp_approx__.Median, tol = abs(.01*maxlp))) |> 
    ungroup()
  #  Find the best ones:
  print("grouping, now binding")
  better = output |> filter(same_as_best == TRUE)
  worse = output  |> filter(same_as_best == FALSE)
  
  saveRDS(better, file = paste0(vb_summary_path, "summary_of_BEST_igt_subgroup_", individuals, "_VBinit_vb.rds"))
  
  saveRDS(worse, file = paste0(vb_summary_path, "summary_of_NOTBEST_igt_subgroup_", individuals, "_VBinit_vb.rds"))
  
  print("saved, better and worse")
  
  
  ################# MERGE THE RESULTS BY MODEL #################
  
  modelsets  = better |> select(model)|> unique()
  # print("MODELSETS:")
  # print(modelsets)
  
  VB_RUNS = NULL
  for(model_grab in 1:nrow(modelsets)){
    # print(paste0("START MODEL_GRAB: ", model_grab))
    f_use      <- better |> 
      filter( model           == modelsets[model_grab,"model"][[1]]) |> pull(filename)
    # print(f_use)
    print("Pulling by models that were 'better'")
    
    combined_samples = NULL
    for(file_index in f_use){ # f_use # f_use[1:min(15, length(f_use))]
      # here is where we discard burn in and concatenate all the MCMC runs.
      print(file_index)
      
      file_in_use = readRDS(paste0(vb_output_path,  file_index ))
      # print(file_in_use)
      
      if(is.null(combined_samples)){# first one
        combined_samples = file_in_use
      }else{# already exists
        for(ind in 1:length(file_in_use$parVals)){
          if(!is.matrix(file_in_use$parVals[[ind]])){# it's a vector
            combined_samples$parVals[[ind]] = c(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]]) 
          }else{#2-dims
            combined_samples$parVals[[ind]] = rbind(combined_samples$parVals[[ind]],file_in_use$parVals[[ind]])
          }
        }# end for over parVals
        
        if (!is.null(file_in_use$fit) && inherits(file_in_use$fit, "stanfit")) {
          if (!is.null(combined_samples$fit) && inherits(combined_samples$fit, "stanfit")) {
            combined_samples$fit = sflist2stanfit(list(combined_samples$fit, file_in_use$fit))
            # print(paste0("combined_samples$fit has values, with size: ", length(combined_samples$fit)))
          } else {
            combined_samples$fit = file_in_use$fit  # First valid assignment
            print("Initializing combined_samples$fit with first valid stanfit object")
          }
        } else {
          print(paste0("Skipping file: ", file_index, " as it does not contain a valid stanfit object"))
        }
        
      }# end else
      
      # combine them all together in this format as well...
      model           = file_in_use$model
      #check if both chains are the same-ish  
      extracted_samples  = rstan::extract(file_in_use$fit,  permuted=FALSE, inc_warmup=TRUE)
      # print(paste0("Extracted Samples: ", extracted_samples))
      if(dim(extracted_samples)[2]==1){
        VB_RUNS        =   bind_rows(VB_RUNS,
                                     data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
      }else{   # eventually allow an arbitrary number of chains, but for now, just be lazy
        VB_RUNS        =   bind_rows(VB_RUNS,
                                     data.frame(extracted_samples[,1,],chain = 1, filename = file_index, model = modelsets[model_grab,"model"][[1]]),
                                     data.frame(extracted_samples[,2,],chain = 2, filename = file_index, model = modelsets[model_grab,"model"][[1]]))
      }
      
    }
    combined_samples |> saveRDS(file = paste0(vb_summary_path, "merged_VB_output_from_igt_subgroup_", individuals, "_",
                                              modelsets[model_grab,"model"][[1]],
                                              ".rds"))
    extract_ic(combined_samples, ic="looic") |> saveRDS(file = paste0(vb_summary_path, "merged_VB_output_from_igt_subgroup_",
                                                                      individuals, "_",
                                                                      modelsets[model_grab,"model"][[1]],
                                                                      "_justlooic.rds"))
  }
  
  saveRDS(VB_RUNS, file = paste0(vb_summary_path, "samples_from_VB_igt_subgroup_", individuals, "_VBinit_vb_estimate.rds"))
  print("samples saved")
}


## Load from VB Summary a bic_results
load_bic_results_function <- function(vb_file_path, individuals_string_with_space) {
  vb_output_path = paste0(vb_file_path, "vb_output/")
  vb_summary_path = paste0(vb_output_path, "summary/")
  individuals_list = unlist(strsplit(individuals_string_with_space, split = " "))
  individuals_string = paste0(individuals_list, collapse = "") 
  
  model_names = c('vpp', 'orl', 'pvl_delta')
  
  ### Read in all model VB outputs ###
  model_vb_results <- list()
  for (m in 1:length(model_names)) {
    file_to_import <- paste0(vb_summary_path, "merged_VB_output_from_igt_subgroup_", individuals_string, "_igt_", model_names[m], ".rds")
    print(paste0("Importing file: ", file_to_import))
    model_vb_results[[paste0(model_names[m], "_vb")]]  <- readRDS(file = file_to_import)
  }
  
  # Create object to store BIC for each individual in each model
  bic_results <- list()
  
  for (model in model_names) {
    # Update name to be including _top15, as this is how it is saved in our data
    model_name <- paste0(model, "_vb")
    
    # Get a parameters list based on the model examined; exclude the subjID
    parameter_names_list <- colnames(model_vb_results[[model_name]]$allIndPars)[2:ncol(model_vb_results[[model_name]]$allIndPars)]
    loglik_df <- as.data.frame.matrix(model_vb_results[[model_name]]$parVals$log_lik)
    loglik_mean <- apply(loglik_df, 2, mean)
    
    # Create BIC storage df:
    bic_df <- data.frame(matrix(ncol = 1, nrow = 0))
    # get which individuals are in the models:
    individuals_in_model <- model_vb_results[[model_name]]$allIndPars$subjID
    
    # Now we cycle through each individual
    for (i in 1:length(individuals_in_model)) { 
      # Create new dataframe to store log-likelihood and bic score
      individual_log_lik <- loglik_mean[paste0('V',i)]
      # use function to calculate each bic score
      individual_bic <- bic_indiv_calc_fn(individual_log_lik, 100, length(parameter_names_list))
      
      bic_df <- rbind(bic_df, individual_bic)
      rownames(bic_df)[i] <- individuals_in_model[i]
      colnames(bic_df) <- c("bic")
    }
    bic_results[[model]]$individual_bic <- bic_df
  }
  return(bic_results)
}


## select BIC best placement for the iteration of VB
place_vb_bic_function <- function(bic_results) {
  vb_best_fit_set <- list(vpp = c(), orl = c(), pvl_delta = c())
  model_individuals <- bic_results$vpp$individual_bic
  list_of_individuals <- as.numeric(rownames(model_individuals))
  
  for (indiv in 1:nrow(model_individuals)) { 
    # compare columns: vpp_BIC, orl_BIC, pvl_delta_BIC for lowest value:
    vpp_bic <- bic_results$vpp$individual_bic[indiv,]
    orl_bic <- bic_results$orl$individual_bic[indiv,]
    pvl_delta_bic <- bic_results$pvl_delta$individual_bic[indiv,]
    
    ## find model with smallest bic:
    bic_min_index <- which.min(c(vpp_bic, orl_bic, pvl_delta_bic))
    ## find model name that has the lowest bic:
    lowest_bic_model_name <- names(bic_results)[bic_min_index]
    
    ## initial_best_fit_set$vpp # has list of all individual parameters that are the best for vpp
    ## initial_best_fit_set$orl # has list of all individual parameters that are the best for orl
    ## initial_best_fit_set$pvl_delta # has list of all individual parameters that are the best for pvl_delta
    vb_best_fit_set[[lowest_bic_model_name]] <- c(vb_best_fit_set[[lowest_bic_model_name]], list_of_individuals[indiv])
  }
  return(vb_best_fit_set)
}


## Pull VB best results fit
# individuals_string_with_space = '1 2 3 4 5 6 7'
## This function should get the best fits from all three models, given the individuals searched for
# input string of individuals you want to run the VB for and get out new placement of them
pull_vb_summary_fit <- function(data_files_path, vb_file_path, individuals_string_with_space){
  # take input of string with spaces, convert to string with - inbetween to find matching file
  individuals_list = unlist(strsplit(individuals_string_with_space, split = " "))
  individuals_string = paste0(individuals_list, collapse = "-") 
  individuals = paste0(individuals_list, collapse = "") 
  # use individuals_string to find corresponding file in results_combined_text_files
  selected_file_name = paste0('results_individuals_', individuals_string, '.txt')
  
  vb_output_path = paste0(vb_file_path, "vb_output/")
  
  vb_file_list <- list.files(vb_output_path, full.names = FALSE)
  # before running, check if the selected_file_name exists in the vb_output folder already, so I'm not running it multiple times
  # check if any matching exists
  exact_match_pattern <- paste0("^igt_", individuals,"_")
  file_already_exists <- grep(exact_match_pattern, vb_file_list, value = TRUE)
  ###########HERE!
  if (length(file_already_exists) > 0) {
    print("VB has already been run for these individuals. Pulling summary fits...")
    # pull summary fits to calculate BIC
  } else {
    print("This combination of individuals has not been run through estimation. Beginning VB estimation now...")
    run_vb_estimate_function(data_files_path, vb_file_path, selected_file_name)
    run_vb_merge_results_function(vb_file_path, individuals) # currently an error; due to no FIT being pulled (dimensions do not match)
  }
  
  # Once I have all the results run, I need to pull the bic scores:
  bic_results <- load_bic_results_function(vb_file_path, individuals_string_with_space)
  # Now we place all individuals where they belong:
  vb_best_fit_set <- place_vb_bic_function(bic_results)
  return(vb_best_fit_set)
}
# Here we will update some dataframe to have previous placement vs new placement; this will have boolean column to see if there is a change in model placement.


## Stabilizing loop to get a final placement after VB estimates
# keeps track of movement
mcmc_vb_difference_loop_function <- function(mcmc_summary_path, data_files_path, vb_file_path, individuals_string_with_space, N_indivs, passable_percent = 0.2) {
  percent_num <- round(N_indivs*passable_percent, 0)
  
  mcmc_placement <- get_best_placement_function(mcmc_summary_path, individuals_string_with_space)
  model_names <- names(mcmc_placement)
  # populate run results with the beginning of the mcmc ones
  run_results_df <- data.frame()
  for (model in model_names) {
    print(model)
    if (length(mcmc_placement[[model]]) > 0){
      for (i in 1:length(mcmc_placement[[model]])){
        print(mcmc_placement[[model]][i])
        person_id <- as.character(mcmc_placement[[model]][i])
        run_results_df[person_id, 'mcmc_run'] <- model 
      }
    }
  }
  # Ensure rownames are treated as numbers and reorder the dataframe
  run_results_df <- run_results_df[order(as.numeric(rownames(run_results_df))), , drop = FALSE]
  
  run = 1
  changes_in_placement_df <- data.frame()
  diff_count = 1 
  diff_sum = N_indivs
  # we will need to update the placement after the vb init:
  placed_vb_update = mcmc_placement
  
  while (diff_sum > percent_num) {
    print("Starting while loop")
    # for the different models of mcmc
    # do the vb runs for the mcmc placement 
    for (model in model_names) {
      individuals_string_with_space = placed_vb_update[[model]]
      ## Need to skip the models that just have 1 person or less in it:
      if (length(individuals_string_with_space) <= 1 ) {
        run_results_df[individuals_string_with_space, paste0('vb_run_', run)] <- model
      }
      else {
        vb_placement <- pull_vb_summary_fit(data_files_path, vb_file_path, individuals_string_with_space)
        for (model in model_names) {
          if (length(vb_placement[[model]])>0){
            for (i in 1:length(vb_placement[[model]])){
              person_id <- as.character(vb_placement[[model]][i])
              run_results_df[person_id, paste0('vb_run_', run)] <- model 
            }
          }
        }
      }
    }
    run = run + 1
    vb_placement_new = list()
    
    for (row in 1:nrow(run_results_df)) {
      indiv_id = as.character(rownames(run_results_df)[row])
      n_columns <- ncol(run_results_df)
      # compare the last two columns of run_results_df
      old_results <- run_results_df[row, n_columns-1] 
      new_results <- run_results_df[row, n_columns]
      
      # Add new placement of individuals:
      vb_placement_new[[new_results]] <- c(vb_placement_new[[new_results]], indiv_id)
      
      # track changes as boolean
      if (new_results == old_results) {
        changes_in_placement_df[indiv_id, paste0('diff_', diff_count)] = 0
      }
      else{
        changes_in_placement_df[indiv_id, paste0('diff_', diff_count)] = 1
      }
    }

    diff_sum <- sum(changes_in_placement_df[paste0('diff_', diff_count)])
    diff_count = diff_count +1
    print(paste0("Run #", run, " has a difference sum of ", diff_sum))
    
    # update old placement with new placement to continue loop
    placed_vb_update = vb_placement_new
  }
  # Return the final placement of individuals
  final_dataframes <- list('placements' = run_results_df, 'difference' = changes_in_placement_df, 'final_placement' = placed_vb_update)
  return(final_dataframes)
}

