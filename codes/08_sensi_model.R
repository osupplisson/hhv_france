#This code performs the MRP and pseudo-BMA approach
# Importing packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")
library(loo)

# Import data ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)

#Otherwise some joint might take foreever
spatial_objects$list_dpt$geometry <- NULL
set_virus <- c("hsv1", "hsv2", "vzv", "ebv", "cmv")
type_sensi <- c('intercept')

inla.setOption("num.threads", '12:2')
inla.setOption(inla.timeout = 60 * 15)

# Fitting loop --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path_to_fit_baseline <- path_to_fit
for (virus in set_virus) {
  #Loeading names of all fitted model
  list_model_baseline <- purrr::map(base::list.files(
    path_to_fit_baseline,
    pattern = ".RDS",
    all.files = T,
    full.names = FALSE
  ),
  function(x) {
    if (str_detect(x, virus)) {
      x
    }
  }) %>%
    discard(is.null) %>%
    # Set names
    purrr::set_names()
  #removing.RDS and _baseline for proper naming
  names(list_model_baseline) <-
    str_remove_all(names(list_model_baseline), ".RDS")
  names(list_model_baseline) <-
    str_remove(names(list_model_baseline), "_baseline")
  
  #Going through all fitted model and fitting the exact same model but replacing the logit link by a probit one
  for (model in 1:length(list_model_baseline)) {
    for (sensitivity_model in type_sensi) {
      print(paste0(model, "/", length(list_model_baseline)))
      #Name of the model
      fit_name <- names(list_model_baseline)[[model]]
      
      
      #Path to save
      path_to_fit <-
        paste0(
          "hhv_france/clean_data/output/sensitivity_analyses/",
          sensitivity_model,
          "/fit/"
        )
      path_save <- paste0(path_to_fit, fit_name, ".RDS")
      
      print(fit_name)
      print(sensitivity_model)
      print(path_save)
      fit_name_old <- fit_name
      fit_name <- paste0(fit_name_old, "_", sensitivity_model)
      #If not fitted: fit otherwise pass
      if (!file.exists(path_save)) {
        #Loading old model
        fit <-
          readRDS(paste0(path_to_fit_baseline, list_model_baseline[[model]]))
        
        if (sensitivity_model == 'probit') {
          fit$.args$control.family[[1]]$control.link$model <- "probit"
        } else if (sensitivity_model == 'intercept') {
          fit$.args$control.fixed$prec <- 0.1
          fit$.args$control.fixed$prec.intercept <- 0.1
        } else if (sensitivity_model == 'probit_and_fixed') {
          fit$.args$control.fixed$prec <- 0.1
          fit$.args$control.fixed$prec.intercept <- 0.1
          fit$.args$control.family[[1]]$control.link$model <-
            "probit"
        }
        results <- NULL
        
        fit_name <- names(list_model_baseline)[[model]]
        
        #Virus name
        virus_name <- str_split(fit_name, pattern = "\\_")[[1]][2]
        #Data from the fit
        space_effect <- str_split(fit_name, pattern = "\\_")[[1]][3]
        print(virus_name)
        print(space_effect)
        if (str_detect(space_effect, "queen")) {
          graph_bym2_district <- inla_graph_districts_queen
        } else if (str_detect(space_effect, "delaunay")) {
          graph_bym2_district <- inla_graph_districts_delaunay
        } else if (str_detect(space_effect, "soi")) {
          graph_bym2_district <- inla_graph_districts_soi
        } else if (str_detect(space_effect, "gabriel")) {
          graph_bym2_district <- inla_graph_districts_gabriel
        } else if (str_detect(space_effect, "relative")) {
          graph_bym2_district <- inla_graph_districts_relative
        } else if (str_detect(space_effect, "nb")) {
          graph_bym2_district <-
            eval(parse(text = paste0(
              "inla_graph_districts_",
              space_effect
            )))
        }
        inla.pardiso.check()
        #Try to fit
        try(results <- inla.rerun(fit))
        
        #Checking results
        #Refit if needed
        if (!is.null(results)) {
          if (function_check(results) == T) {
            rerun <- 1
            to_run <- T
            max_run <- 6
            while (to_run == T & rerun <= max_run) {
              if (rerun > 3) {
                print("Trying with GSL optimiser")
                results$.args$control.inla$optimiser <- "gsl"
                results$.args$control.inla$tolerance <- 0.0000001
              }
              print("Probably some issues in the fit, rerun")
              results <- try(inla.rerun(results))
              print(results)
              if (class(results) != "try-error") {
                print("Refit:sucess")
                ## If refit worked: either it's ok and stop or we try to refit if rerun<max_rerun
                if (function_check(fit) == T) {
                  print("Still issue, run again")
                  # If issue=> still refit
                  to_run <- T
                } else {
                  print("Everything is fine")
                  # If ok: STOP RERUN
                  to_run <- F
                }
                # Incrementing rerun
                rerun <- rerun + 1
              } else {
                print("Refit:Failure")
                to_run <- F
                results <- NULL
                print("FAILURE :(")
              }
              
              if (rerun > max_run) {
                print("Number of max rerun reached")
                to_run <- F
                results <- NULL
                print("FAILURE :(")
              }
            }
          } else {
            print("INLA did its job :-)")
          }
        } else {
          # If initial fit NA => Failure
          fit <- NULL
          print("FAILURE :(")
        }
        
        
        if (class(results) == "inla") {
          results$data_fit <- fit$data_fit
          results$fomula_good_format <- fit$fomula_good_format
          print("Computing CV")
          results$loocv <-
            inla.group.cv(results,
                          num.level.sets = -1,
                          strategy = "posterior")
          for (x in c(15, 25)) {
            print(x)
            namex <- paste0("lgocv.m", x)
            results[[namex]] <-
              inla.group.cv(
                results,
                num.level.sets = x,
                strategy = "posterior",
                size.max = x
              )
          }
        } else{
          results <- NULL
        }
        # Save the fit
        saveRDS(results,
                path_save)
        rm("results")
      } else{
        print("Model already available")
      }
    }
  }
}


#BMA ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/output/results_bma.rda")
rm(
  "hsv1_bma_analysis",
  "hsv2_bma_analysis",
  "vzv_bma_analysis",
  "cmv_bma_analysis",
  "ebv_bma_analysis"
)


# Applying the function to all viruses --------------------------------------------------------------------------------------------------------------------------------------------------------------
for (virus.for.loop in set_virus) {
  for (sensi in type_sensi) {
    name.for.loop <-
      paste0(virus.for.loop, "_bma_analysis_", sensi)
    path_to_fit <- paste0("hhv_france/clean_data/output/sensitivity_analyses/", sensi, "/fit/")
    
    output_saved <- base::list.files(
      str_replace(path_to_fit, "fit", "post_fit"),
      pattern = ".RDS",
      all.files = T,
      full.names = FALSE
    )
    
    print(path_to_fit)
    # If not available, perform BMA
    if (!paste0(name.for.loop, ".RDS") %in% output_saved) {
      print("Perform BMA for:")
      print(virus.for.loop)
      bma_draw <- summary_analyses(virus.loop = virus.for.loop,
                                   n_draws = 10000)
      
      
      print(bma_draw$bma_results_summary)
      
      # Save
      saveRDS(bma_draw,
              paste0(
                str_replace(path_to_fit, "fit", "post_fit"),
                name.for.loop,
                ".RDS"
              ))
      
      rm("bma_draw")
      gc(full = T)
    }
  }
}
