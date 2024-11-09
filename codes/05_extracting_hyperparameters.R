#This code extract model hyperparameters
# Importing packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")

# Load data ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)

load("hhv_france/clean_data/output/results_bma.rda")


# Extracting hyperparameters values ---------------------------------------------------------------------------------------------------------------------------------------------------------------
exporting_hyperparameters <- function(set_models) {
  #Vector to list
  list_model <- purrr::map(
    set_models$model,
    function(x) {
      x
    }
  ) %>%
    # Set names
    purrr::set_names()
  
  
  # Import
  list_model <- list_model %>%
    purrr::map(
      function(x) {
        readRDS(str_replace_all(paste0(path_to_fit, x), "/", "//"))
      }
    ) 
  #Extracting parameters
  extracting_hyper_param <- purrr::map(
    names(list_model),
    function(model_name) {
      tmp <- list_model[[model_name]]

      summary_hyper <- tmp$summary.hyperpar %>%
        as.data.frame() %>%
        rownames_to_column(var = "param") %>%
        mutate(model = model_name) %>%
        extraction_model_name()

      summary_intercept <- tmp$summary.fixed %>%
        as.data.frame() %>%
        rownames_to_column(var = "param") %>%
        mutate(model = model_name) %>%
        extraction_model_name()

      plyr::rbind.fill(
        summary_hyper,
        summary_intercept
      ) %>%
        virus_extraction(., "model")
    }
  ) %>%
    # Binding everything
    do.call(rbind, .)
  return(extracting_hyper_param)
}



# Exporting hyperramaters -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Applying the function to all viruses --------------------------------------------------------------------------------------------------------------------------------------------------------------
output_saved <- base::list.files(
  "hhv_france/clean_data/output/post_fit_analyses/",
  pattern = ".rda",
  all.files = T,
  full.names = FALSE
)

if(! "hyperparameters.rda" %in% output_saved){
  
hyperparameters_hsv1 <- exporting_hyperparameters(hsv1_bma_analysis$list_model_to_import)
hyperparameters_hsv2 <- exporting_hyperparameters(hsv2_bma_analysis$list_model_to_import)
hyperparameters_vzv <- exporting_hyperparameters(vzv_bma_analysis$list_model_to_import)
hyperparameters_ebv <- exporting_hyperparameters(ebv_bma_analysis$list_model_to_import)
hyperparameters_cmv <- exporting_hyperparameters(cmv_bma_analysis$list_model_to_import)

# Saving ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gc()
rm(list=setdiff(ls(), paste0("hyperparameters_",c("hsv1", "hsv2", "cmv", "vzv", "ebv"))))
save.image("hhv_france/clean_data/output/post_fit_analyses/hyperparameters.rda")
}
