#This code performs the MRP and BMA approach
# Importing packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")
library(loo)

# Import data ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)

#Otherwise some joint might take foreever
spatial_objects$list_dpt$geometry <- NULL



# Calling the functions --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SEE: https://discourse.mc-stan.org/t/implementing-model-averaging/6082/4



list_strata <- list(
  # Whole France
  c("all"),
  c("sex_numeric"),
  c("age_numeric"),
  c("year_numeric"),
  
  # By districts, region, area
  c("district_numeric"),
  c("nom_region"),
  c("france_area"),
  
  # By area
  c("sex_numeric", "france_area"),
  c("age_numeric", "france_area"),
  c("year_numeric", "france_area"),
  
  # By region
  c("sex_numeric", "nom_region"),
  c("age_numeric", "nom_region"),
  c("year_numeric", "nom_region"),
  
  # By district
  c("sex_numeric", "district_numeric"),
  c("age_numeric", "district_numeric"),
  c("year_numeric", "district_numeric"),
  
  
  # Sex x Age
  c("sex_numeric", "age_numeric"),
  c("sex_numeric", "age_numeric", "france_area"),
  c("sex_numeric", "age_numeric", "district_numeric"),
  c("sex_numeric", "age_numeric", "nom_region"),
  
  
  # Sex x Year
  c("sex_numeric", "year_numeric"),
  c("sex_numeric", "year_numeric", "france_area"),
  c("sex_numeric", "year_numeric", "district_numeric"),
  c("sex_numeric", "year_numeric", "nom_region"),
  
  
  # Comparing with https://www.cambridge.org/core/journals/epidemiology-and-infection/article/seroprevalence-of-cytomegalovirus-infection-in-france-in-2010/FBEF04F5E2A1F155B5AD515D67D85070
  c("age_between_15_and_49_yo"),
  c("age_between_15_and_49_yo", "year_numeric", "france_area"),
  c("age_between_15_and_49_yo", "district_numeric"),
  c("age_between_15_and_49_yo", "nom_region"),
  
  
  
  c("age_class_for_comparisions_with_antona"),
  c(
    "age_class_for_comparisions_with_antona",
    "year_numeric",
    "france_area"
  ),
  c(
    "age_class_for_comparisions_with_antona",
    "district_numeric"
  ),
  c("age_class_for_comparisions_with_antona",
    "nom_region"),
  
  
  c("france_area", "age_numeric"),
  c("france_area", "sex_numeric"),
  c("france_area", "sex_numeric", "age_numeric"),
  
  c("age_class_yousuf"),
  c("age_class_yousuf", "france_area"),
  c("age_class_yousuf", "year_numeric", "france_area")
)


## Vector used for dynamically select cols
c_one_of <- c(
  "all",
  "sex_numeric",
  "age_numeric",
  "year_numeric",
  "district_numeric",
  "nom_region",
  "france_area",
  "age_between_15_and_49_yo",
  "age_class_for_comparisions_with_antona",
  "age_class_yousuf"
)



summary_analyses <- function(virus.loop = "hsv1",
                             n_draws = 5,
                             type_log_score = "lgocv.m15",
                             specific_path = NULL,
                             input_list_weights = NULL,
                             using_weighted_bma_input = FALSE) {
  print("Importing models")
  # Import all models -----------------------------------------------------------------------------------------------------------------------------------------------------------
  ## Create a list with all models that should be considered for averaging
  # @name_model: type of model to consider
  # @virus: virus to consider
  # @package: package to consider
  ## First restrict model list
  ## Then name the list
  ## Then import the model from each list element
  ## Output: A list with name *model name* and list element corresponding fit
  list_model <- purrr::map(base::list.files(
    path_to_fit,
    pattern = ".RDS",
    all.files = T,
    full.names = FALSE
  ),
  function(x) {
    if (str_detect(x, virus.loop)) {
      x
    }
  }) %>%
    # Set names
    purrr::set_names() %>%
    purrr::discard(is.null)
  
  print(list_model)
  
  # Import
  list_model <- list_model %>%
    purrr::map(function(x) {
      readRDS(str_replace_all(paste0(path_to_fit, x), "/", "//"))
    }) %>%
    purrr::discard(is.null)
  
  
  print(list_model)
  if (str_detect(path_to_fit, "sensitivity")) {
    print("Removing models with issues")
    list_model <- list_model %>%
      imap(function(.x, .y) {
        if (length(.x) == 1) {
          NULL
        } else{
          .x
        }
      }) %>%
      purrr::discard(is.null)
  }
  
  
  print("Weights of each stratum")
  #All models have the same data, we just pick from the first model
  observed_strata <- list_model[[1]]$data_fit %>%
    filter(name == "observed") %>%
    dplyr::select(sex_numeric,
                  age_numeric,
                  year_numeric,
                  district_numeric,
                  status = name,
                  id_data)
  
  strata_weights <- list_model[[1]]$data_fit %>%
    filter(name == "predicted") %>%
    dplyr::select((-id_data)) %>%
    left_join(
      observed_strata,
      by = c(
        "sex_numeric",
        "age_numeric",
        "year_numeric",
        "district_numeric"
      )
    ) %>%
    mutate(population_weights = N / sum(N)) %>%
    dplyr::select(population_weights,
                  id_data) %>%
    filter(!is.na(id_data)) %>%
    arrange(id_data)
  
  #Share of stratum with at least 1 observation:
  sum(strata_weights$population_weights)
  
  print("Assessment")
  
  assessment <- list_model %>%
    imap(function(.x, .y) {
      id_obs <- .x$data_fit %>%  filter(name == "observed") %>% .$id_data
      
      if (!all.equal(id_obs, strata_weights$id_data)) {
        print("Issue in the index, STOP")
        stop()
      }
      #Final table
      tmp <- tibble(
        model = .y,
        dic = .x$dic$dic,
        waic = .x$waic$waic
      )
      #Adding metrics
      for (x in c("loocv", paste0("lgocv.m", c(15, 25)))) {
        pd <- .x[[x]]$cv[id_obs]
        ls <- -sum(log(pd))
        weighted_log_pd <-
          pd * log(strata_weights$population_weights)
        weighted_ls <- -sum(weighted_log_pd)
        tmp[[paste0("ls_", x)]] <- ls
        tmp[[paste0("weighted_ls_", x)]] <- weighted_ls
      }
      
      tmp
    }) %>%
    do.call(rbind, .) %>%
    dplyr::select(model,
                  dic,
                  waic,
                  starts_with("ls"),
                  starts_with("weighted_ls")) %>%
    pivot_longer(cols = c(-model)) %>%
    group_by(name) %>%
    arrange(value) %>%
    mutate(rank = row_number())  %>%
    pivot_wider(names_from = name,
                values_from = c(value, rank)) %>%
    virus_extraction(., "model") %>%
    extraction_model_name() %>%
    rename_all(function(.) {
      str_remove(., "value\\_")
    })
  
  print(assessment)
  print(assessment$rank_weighted_ls_lgocv.m15)
  print(assessment$rank_ls_lgocv.m15)
  
  print("Stacking Weights")
  
  function_for_weights <- function(metric = "lgocv.m15",
                                   input_model_list = list_model,
                                   using_weighted_bma = using_weighted_bma_input) {
    # Weights associated with each model
    log_gcpo <- purrr::map(names(input_model_list),
                           function(model.loop) {
                             # Extract results
                             model_results <-
                               input_model_list[[model.loop]]
                             
                             id_obs <-
                               model_results$data_fit %>%
                               filter(name == "observed") %>%
                               .$id_data
                             
                             # Extract Log-metric
                             vector <- model_results[[metric]]$cv
                             vector <- vector[id_obs]
                             log_gcpo <- log(vector)
                             
                             log_gcpo <-
                               tibble(log_gcpo = log_gcpo,
                                      id_data = id_obs) %>%
                               mutate("model" = model.loop,
                                      id = row_number())
                             
                             # Adding survey weights
                             if (using_weighted_bma == T) {
                               print("Weighting by population weigths")
                               
                               #Keeping only observed strata
                               log_gcpo <- log_gcpo %>%
                                 left_join(strata_weights,
                                           by = "id_data")
                               
                               #Weigthing gcpo/cpo
                               log_gcpo <- log_gcpo %>%
                                 mutate(log_gcpo = log_gcpo * population_weights)
                             }
                             log_gcpo <- log_gcpo %>%
                               dplyr::select(log_gcpo,
                                             id,
                                             model)
                           }) %>%
      # Binding all log_lgocv
      do.call(rbind, .) %>%
      # Pivoting
      pivot_wider(names_from = model,
                  values_from = log_gcpo) %>%
      # Removing ID
      dplyr::select(-id)
    
    # Feed the matrix of lgocv to the loo-package
    weights <- tibble(colnames(log_gcpo),
                      stacking_weights(as.matrix(log_gcpo)))
    
    # Computing loo score
    N <- nrow(as.matrix(log_gcpo))
    K <- ncol(as.matrix(log_gcpo))
    negative_log_score_loo <- function(w) {
      sum <- 0
      for (i in 1:N) {
        sum <- sum + log(exp(as.matrix(log_gcpo)[i,]) %*% w)
      }
      # We look for minimum neg score
      return(-as.numeric(sum))
    }
    
    
    
    # Adding names
    colnames(weights) <- c("model",
                           "stacking")
    
    stacking_score <- negative_log_score_loo(w = weights$stacking)
    #Adding stacking score
    weights$stacking_score <- stacking_score
    
    
    return(weights)
  }
  
  function_for_filtering_out_model <- function(metric = "loocv",
                                               input_model_list = list_model) {
    print(metric)
    weights <- function_for_weights(metric = metric,
                                    input_model_list = input_model_list)
    model_to_drop_loop <- weights %>%
      arrange(desc(stacking)) %>%
      mutate(cum = cumsum(stacking)) %>%
      filter(cum > 0.99) %>%
      filter(stacking < 0.01)
    
    if (nrow(model_to_drop_loop) == 0) {
      model_dropped <- model_to_drop_loop
    }
    # Filtering out models up until BMA is best or selecting only one model
    print("Filtering out model with negligible weights")
    while (nrow(model_to_drop_loop) > 0) {
      print('Removing the following model:')
      print(model_to_drop_loop)
      if (!exists("model_dropped")) {
        model_dropped <- model_to_drop_loop
      } else{
        model_dropped <- rbind(model_dropped,
                               model_to_drop_loop)
      }
      
      print('Models dropped so far:')
      print(model_dropped)
      
      #Removing model
      for (m in model_to_drop_loop$model) {
        input_model_list[[m]] <- NULL
      }
      
      rm("model_to_drop_loop")
      if (length(input_model_list) > 1) {
        #Update weights
        weights <- function_for_weights(
          metric = metric,
          input_model_list = input_model_list,
          using_weighted_bma = using_weighted_bma_input
        )
        #Do we have models to drop ?
        model_to_drop_loop <- weights %>%
          arrange(desc(stacking)) %>%
          mutate(cum = cumsum(stacking)) %>%
          filter(cum > 0.99) %>%
          filter(stacking < 0.01)
      } else{
        weights <- tibble(
          model = names(input_model_list),
          stacking = 1,
          stacking_score = assessment %>% filter(model == names(input_model_list)) %>% .$"lgocv.m15"
        )
        model_to_drop_loop <- tibble()
      }
    }
    
    return(list("weights" = weights,
                "model_dropped" = model_dropped))
  }
  
  if (is.null(input_list_weights)) {
    list_weights <- list()
    for (z in c("loocv", paste0("lgocv.m", c(15, 25)))) {
      list_weights[[z]] <- function_for_filtering_out_model(metric = z,
                                                            input_model_list = list_model)
    }
    rm("z")
  } else{
    print('Using previously computed lists of weights to gain time')
    print("Using weights:")
    print(type_log_score)
    list_weights <- input_list_weights
  }
  
  
  if (is.null(specific_path)) {
    print("Chosen weights")
    print(type_log_score)
    weights <- list_weights[[type_log_score]]$weights
    model_dropped <- list_weights[[type_log_score]]$model_dropped
    print("Remaining models")
    print(weights)
    if (!is.null(model_dropped)) {
      print("Dropping models")
      for (m in model_dropped$model) {
        list_model[[m]] <- NULL
      }
    } else{
      "No model dropped"
    }
    
    
    weights <- weights %>%
      mutate(type_log_score = type_log_score)
    
    # Extracting virus and adding number of draws to keep for each model
    # ROUNDING WEIGHTS HAPPENS AT THIS STEP !
    weights <- weights %>%
      virus_extraction(., "model") %>%
      mutate(n_to_draw =  as.integer(round(stacking * n_draws))) %>%
      extraction_model_name()
  } else{
    ##If instead of BMA we wantto use th best model
    if (specific_path == "sensitivity_best_model") {
      print("Watch out!")
      print("Specific 'sensitivity' path")
      print("BMA performed only on the best model according to log-lgocv score")
      weights <- assessment %>%
        mutate(n_to_draw = ifelse(rank_ls_lgocv.m15 == 1,
                                  n_draws,
                                  0)) %>%
        dplyr::select(model,
                      virus,
                      starts_with("label"),
                      n_to_draw) %>%
        arrange(desc(n_draws))
      
      print(weights)
      
      
      print("Dropping models")
      for (m in names(list_model)) {
        if (m != weights %>% filter(n_to_draw == n_draws) %>% .$model) {
          list_model[[m]] <- NULL
        }
      }
    }
    else{
      ##If instead of BMA we wantto use th best model
      if (specific_path == "sensitivity_best_weighted_model") {
        print("Watch out!")
        print("Specific 'sensitivity' path")
        print("BMA performed only on the best model according to weighted log-lgocv score")
        weights <- assessment %>%
          mutate(n_to_draw = ifelse(rank_weighted_ls_lgocv.m15 == 1,
                                    n_draws,
                                    0)) %>%
          dplyr::select(model,
                        virus,
                        starts_with("label"),
                        n_to_draw) %>%
          arrange(desc(n_draws))
        
        print(weights)
        
        
        print("Dropping models")
        for (m in names(list_model)) {
          if (m != weights %>% filter(n_to_draw == n_draws) %>% .$model) {
            list_model[[m]] <- NULL
          }
        }
      }
    }
  }
  
  # Saving list of models with weights > 0 to import for PP check and hyperparameters
  list_model_to_stack <- weights %>%
    filter(n_to_draw > 0L)
  
  print(list_model_to_stack)
  
  print("Filtering out models with 0 draws")
  # Filtering out from the list models with n_to_draw <= 0
  list_model <-
    list_model[names(list_model) %in% list_model_to_stack$model]
  
  ##### weights contains the final weights/model
  print("Draws for PS and BMA")
  # Draws for post-stratification -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## For each fitted model:
  ### 1- Extract the fit
  ### 2- Extract data corresponding to the poststratification using the true population of each stratum and keep their ids
  ### 3- Draws the mu for each stratum proportionally to the weights
  ### 5- Stack everything
  all_draws <- list_model %>%
    imap(function(.x, .y) {
      # Extract the fit
      tmp <- .x
      
      # Extract the data used for the fit, keep only rows corresponding to the strata
      strata <- tmp$data_fit %>%
        dplyr::select(-outcome,
                      -one_of("id_data")) %>%
        mutate(id_row = row_number()) %>%
        filter(name == "predicted") %>%
        mutate(strata_id = paste0(sex_numeric,
                                  age_numeric,
                                  year_numeric,
                                  district_numeric))
      
      # Extract the number of draws to perform for each model
      n_to_draw <- weights %>%
        filter(model == .y) %>%
        .$n_to_draw
      
      
      # They are all >0
      # Draws from the joint posterior
      draw_approx_joint <-
        inla.posterior.sample(n = n_to_draw,
                              result = tmp,
                              seed = 123)
      
      # Function to extract the parameter of the binomial
      sim.p <- function() {
        plogis(Predictor)
      }
      if (str_detect(.y, "probit")) {
        sim.p <- function() {
          pnorm(Predictor)
        }
      }
      print(sim.p)
      
      
      # Extract the parameter of the binomial
      mu_draws <-
        inla.posterior.sample.eval(sim.p,
                                   draw_approx_joint)
      
      # Keeping only rows from the strata in the fit
      mu_draws <-
        mu_draws[strata$id_row, ] %>%
        as_tibble()
      if (n_to_draw == 1) {
        # If n_to_draw=1, changing the colname for consistency with n_to_draw>1
        colnames(mu_draws) <- "sample:1"
      }
      
      # Linking draws and strata information
      data_draws <- mu_draws %>%
        cbind(strata) %>%
        pivot_longer(
          cols = starts_with("sample:"),
          names_to = "sample",
          values_to = "p"
        ) %>%
        mutate(model = .y) %>%
        virus_extraction(., "model") %>%
        dplyr::select(
          virus,
          model,
          strata_id,
          sex_numeric,
          age_numeric,
          district_numeric,
          year_numeric,
          N,
          sample,
          p
        ) %>%
        left_join(
          spatial_objects$list_dpt %>% dplyr::select(district_numeric, nom_region, insee_reg),
          by = "district_numeric"
        ) %>%
        mutate(all = "all") %>%
        function_france_area() %>%
        left_join(age_vector %>% dplyr::select(age_class, age_numeric), by = "age_numeric") %>%
        # Adding age categories for comparisions with https://www.cambridge.org/core/journals/epidemiology-and-infection/article/seroprevalence-of-cytomegalovirus-infection-in-france-in-2010/FBEF04F5E2A1F155B5AD515D67D85070
        mutate(
          age_between_15_and_49_yo = ifelse(
            age_class %in% c(
              "(14,19]",
              "(19,24]",
              "(24,29]",
              "(29,34]",
              "(34,39]",
              "(39,44]",
              "(44,49]"
            ),
            "[15-49]",
            "<15 or >49 yo"
          ),
          age_class_for_comparisions_with_antona = case_when(
            age_class %in% c("(14,19]", "(19,24]") ~ "(14,24]",
            age_class %in% c("(24,29]", "(29,34]") ~ "(25,34]",
            age_class %in% c("(34,39]", "(39,44]", "(44,49]") ~ "(34,49]",
            TRUE ~ age_class
          ),
          age_class_yousuf = case_when(
            age_class %in% c("[0,4]", "(4,9]" , "(9,14]", "(14,19]") ~ "<20",
            age_class %in% c("(19,24]", "(24,29]") ~ "[20,30)",
            age_class %in% c("(29,34]", "(34,39]") ~ "[30,40)",
            age_class %in% c("(39,44]", "(44,49]") ~ "[40,50)",
            TRUE ~ ">49"
          )
        ) %>%
        dplyr::select(-age_class)
    }) %>%
    # Binding everything
    do.call(rbind, .)
  
  print("PS")
  # Perform the post-stratification -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # @df_draws: output from function_to_compute_weights
  # @group_to_stratify: strata
  post_stratification_value <- purrr::map(list_strata,
                                          function(group_to_stratify) {
                                            group_to_stratify
                                          }) %>%
    set_names() %>%
    purrr::map(function(group_to_stratify) {
      all_draws %>%
        # For each model, draws, and subpopulation
        group_by(model,
                 sample,
                 !!!rlang::syms(group_to_stratify)) %>%
        summarise(# Computing the total population
          N_tot = sum(N),
          # Post-stratifying
          p_ps = 100 * sum(N * p / sum(N))) %>%
        ungroup() %>%
        # Computing the initial number of draws/model
        group_by(model,
                 !!!rlang::syms(group_to_stratify)) %>%
        mutate(n_draws = n()) %>%
        ungroup() %>%
        # Adding N_tot
        virus_extraction(., "model") %>%
        dplyr::select(virus,
                      model,
                      n_draws,
                      sample,
                      one_of("model",
                             c_one_of),
                      N_tot,
                      p_ps)
    })
  
  ## Adding some data frame to compute differences b/w strata
  ####Global difference
  post_stratification_value$diff_sex <-
    post_stratification_value$sex_numeric
  post_stratification_value$diff_age <-
    post_stratification_value$age_numeric
  post_stratification_value$diff_year <-
    post_stratification_value$year_numeric
  
  ####Difference by area
  post_stratification_value$diff_sex_france <-
    post_stratification_value$`c("sex_numeric", "france_area")`
  post_stratification_value$diff_age_france <-
    post_stratification_value$`c("age_numeric", "france_area")`
  post_stratification_value$diff_year_france <-
    post_stratification_value$`c("year_numeric", "france_area")`
  
  
  # post_stratification_value$special_district <-
  #   post_stratification_value$district_numeric
  # post_stratification_value$special_region <-
  #   post_stratification_value$nom_region
  
  post_stratification_value$pairwise_age <-
    post_stratification_value$age_numeric
  post_stratification_value$pairwise_year <-
    post_stratification_value$year_numeric
  post_stratification_value$pairwise_district <-
    post_stratification_value$district_numeric
  post_stratification_value$pairwise_region <-
    post_stratification_value$nom_region
  post_stratification_value$pairwise_area <-
    post_stratification_value$france_area
  
  post_stratification_value$pairwise_age_france <-
    post_stratification_value$`c("age_numeric", "france_area")`
  post_stratification_value$pairwise_year_france <-
    post_stratification_value$`c("year_numeric", "france_area")`
  
  post_stratification_value$pairwise_age_sex <-
    post_stratification_value$`c("sex_numeric", "age_numeric")`
  post_stratification_value$pairwise_age_sex_france <-
    post_stratification_value$`c("sex_numeric", "age_numeric", "france_area")`
  
  post_stratification_value$pairwise_sex_age <-
    post_stratification_value$`c("sex_numeric", "age_numeric")`
  post_stratification_value$pairwise_sex_age_france <-
    post_stratification_value$`c("sex_numeric", "age_numeric", "france_area")`
  
  
  post_stratification_value$pairwise_france_age <-
    post_stratification_value$`c("france_area", "age_numeric")`
  post_stratification_value$pairwise_france_sex <-
    post_stratification_value$`c("france_area", "sex_numeric")`
  post_stratification_value$pairwise_france_sex_age <-
    post_stratification_value$`c("france_area", "sex_numeric", "age_numeric")`
  
  numero_paris <-
    spatial_objects$list_dpt %>%
    filter(nom == "Paris") %>%
    .$district_numeric
  region_paris <-
    spatial_objects$list_dpt %>%
    filter(nom == "Paris") %>%
    .$nom_region
  
  print("Summary BMA")
  bma_results_summary <- post_stratification_value %>%
    purrr::imap(# If not difference: compute classical seroprevalence
      ~ if (!str_detect(.y, "diff") &
            !str_detect(.y, "special") &
            !str_detect(.y, "pairwise")) {
        .x %>%
          mutate(value = p_ps) %>%
          group_by_at(vars(one_of(c_one_of))) %>%
          function_for_summary_draws() %>%
          mutate(type = "expected_seroprev")
      }
      else if (str_detect(.y, "diff")) {
        # Else compute difference and ratio in seroprevalence with previous age class or year
        #First variable to group with for ordering
        var_to_group_by_1 <-
          case_when(str_detect(.y, "france") ~ "france_area",
                    TRUE ~ NA)
        if (is.na(var_to_group_by_1)) {
          var_to_group_by_1 <- NULL
        }
        #Second variable to group with for summary
        var_to_group_by_2 <- case_when(
          str_detect(.y, "sex") ~ "sex_numeric",
          str_detect(.y, "age") ~ "age_numeric",
          str_detect(.y, "year") ~ "year_numeric"
        )
        
        if (!is.null("var_to_group_by_1")) {
          var_to_group_by_2 <- c(var_to_group_by_1,
                                 var_to_group_by_2)
        }
        
        .x %>%
          #/!\  Difference should be computed for a given virus/model/draws
          group_by_at(c("virus", "model", "sample", var_to_group_by_1)) %>%
          arrange_at(var_to_group_by_2) %>%
          mutate(value = p_ps - lag(p_ps)) %>%
          ungroup() %>%
          group_by_at(c("virus", var_to_group_by_2)) %>%
          function_for_summary_draws()
      }
      else if (str_detect(.y, "special")) {
        # Else compute difference in seroprevalence wrt PARIS
        var_to_group_by <- case_when(
          str_detect(.y, "district") ~ "district_numeric",
          str_detect(.y, "region") ~ "nom_region"
        )
        if (str_detect(.y, "district")) {
          ref_df <- .x[.x$district_numeric == numero_paris, ]
          ref_df <-
            ref_df %>% rename(p_ps_ref = p_ps, district_ref = district_numeric)
        } else {
          ref_df <- .x[.x$nom_region == region_paris, ]
          ref_df <-
            ref_df %>% rename(p_ps_ref = p_ps, nom_region_ref = nom_region)
        }
        
        .x %>%
          left_join(ref_df, by = c("virus", "model", "sample")) %>%
          group_by_at(c("virus", "model", "sample")) %>%
          arrange_at(var_to_group_by) %>%
          mutate(value = p_ps - p_ps_ref) %>%
          ungroup() %>%
          group_by_at(c("virus", var_to_group_by)) %>%
          function_for_summary_draws() %>%
          mutate(type = "first_diff")
      }
      else if (str_detect(.y, "pairwise")) {
        #Second variable to group with for summary
        var_to_group_by <- case_when(
          .y %in% c(
            'pairwise_age',
            'pairwise_age_france',
            'pairwise_age_sex',
            'pairwise_age_sex_france'
          ) ~ "age_numeric",
          .y %in% c('pairwise_year',
                    'pairwise_year_france') ~ 'year_numeric',
          .y == 'pairwise_district' ~ "district_numeric",
          .y == 'pairwise_region' ~ "nom_region",
          .y == 'pairwise_area' ~ "france_area",
          .y %in% c('pairwise_sex_age',
                    'pairwise_sex_age_france') ~ 'sex_numeric',
          .y %in% c(
            "pairwise_france_age",
            "pairwise_france_sex",
            "pairwise_france_sex_age"
          ) ~ 'france_area',
          TRUE ~ "NA"
        )
        
        
        .x[["var_to_diff"]] <- .x[[var_to_group_by]]
        
        
        var_to_group_by_2 <- case_when(
          str_detect(.y, "_france") &
            !str_detect(.y, "pairwise_france") ~ 'france_area',
          .y == 'pairwise_age_sex' ~ 'sex_numeric',
          .y == 'pairwise_sex_age' ~ 'age_numeric',
          .y == 'pairwise_france_age' ~ 'age_numeric',
          .y == 'pairwise_france_sex' ~ 'sex_numeric',
          .y == 'pairwise_france_sex_age' ~ 'sex_numeric',
          TRUE ~ "NA"
        )
        
        
        var_to_group_by_3 <- case_when(
          .y == 'pairwise_age_sex_france' ~ 'sex_numeric',
          .y == 'pairwise_sex_age_france' ~ 'age_numeric',
          .y == 'pairwise_france_sex_age' ~ 'age_numeric',
          TRUE ~ "NA"
        )
        
        
        
        .x %>%
          ungroup() %>%
          #/!\  Difference should be computed for a given virus/model/draws
          group_by_at(vars(one_of(
            c(
              "virus",
              "model",
              "sample",
              var_to_group_by_2,
              var_to_group_by_3
            )
          )))  %>%
          summarise(
            value = combn(p_ps, 2, diff),
            var_to_diff = combn(var_to_diff, 2, paste0, collapse = '-')
          ) %>%
          ungroup() %>%
          group_by_at(vars(one_of(
            c(
              "virus",
              "var_to_diff",
              var_to_group_by_2,
              var_to_group_by_3
            )
          ))) %>%
          function_for_summary_draws()
        
      }) %>%
    # Adding variable name
    purrr::imap(~ if ("sex_numeric" %in% colnames(.x)) {
      left_join(.x, sex_vector, by = "sex_numeric")
    } else {
      .x
    }) %>%
    purrr::imap(~ if ("age_numeric" %in% colnames(.x)) {
      left_join(.x, age_vector, by = "age_numeric")
    } else {
      .x
    }) %>%
    purrr::imap(~ if ("year_numeric" %in% colnames(.x)) {
      left_join(.x, year_vector, by = "year_numeric")
    } else {
      .x
    }) %>%
    purrr::imap(~ if ("district_numeric" %in% colnames(.x)) {
      left_join(.x, district_vector, by = "district_numeric")
    } else {
      .x
    })
  
  
  
  
  output <- list(
    "list_model_to_import" = list_model_to_stack,
    "assessment" = assessment,
    #List weights
    "list_weights" = list_weights,
    #Chosen weight + samples
    "weights" = weights,
    "bma_results_summary" = bma_results_summary
  )
  return(output)
}



# Applying the function to all viruses --------------------------------------------------------------------------------------------------------------------------------------------------------------
n_draws_choice <- 10000
for (virus.for.loop in c("hsv1", "hsv2", "vzv", "ebv", "cmv")) {
  name.for.loop <- paste0(virus.for.loop, "_bma_analysis")
  saving_path <-
    paste0("hhv_france/clean_data/output/post_fit_analyses/",
           name.for.loop,
           ".RDS")
  
  if (file.exists(saving_path) == FALSE) {
    print("Perform baseline BMA for:")
    print(virus.for.loop)
    bma_draw <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      type_log_score = "lgocv.m15"
    )
    
    saveRDS(bma_draw, saving_path)
    
    print(bma_draw$bma_results_summary)
  } else{
    print("Already available")
    bma_draw <- readRDS(saving_path)
  }
  
  
  name.for.loop <-
    paste0(virus.for.loop, "_mrp_analysis_best_model_only")
  saving_path <- paste0(
    "hhv_france/clean_data/output/sensitivity_analyses/mrp_best_model_only/",
    name.for.loop,
    ".RDS"
  )
  if (file.exists(saving_path) == FALSE) {
    print("Perform BMA on best model")
    x <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      specific_path = "sensitivity_best_model",
      input_list_weights = bma_draw$list_weights
    )
    saveRDS(x, saving_path)
    
    rm("x")
  } else{
    print("Already available")
  }
  
  
  name.for.loop <-
    paste0(virus.for.loop, "_mrp_analysis_best_weighted_model_only")
  saving_path <- paste0(
    "hhv_france/clean_data/output/sensitivity_analyses/mrp_best_model_only/",
    name.for.loop,
    ".RDS"
  )
  if (file.exists(saving_path) == FALSE) {
    print("Perform BMA on best WEIGHTED model")
    x <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      specific_path = "sensitivity_best_weighted_model",
      input_list_weights = bma_draw$list_weights
    )
    saveRDS(x, saving_path)
    
    rm("x")
  } else{
    print("Already available")
  }
  
  
  
  name.for.loop <- paste0(virus.for.loop, "_wbma_analysis")
  saving_path <- paste0(
    "hhv_france/clean_data/output/sensitivity_analyses/weighted_bma/",
    name.for.loop,
    ".RDS"
  )
  if (file.exists(saving_path) == FALSE) {
    print("Perform BMA with strata weighted by their proportion in the population")
    
    x <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      type_log_score = "lgocv.m15",
      using_weighted_bma_input = TRUE
    )
    saveRDS(x, saving_path)
    rm("x")
  } else{
    print("Already available")
  }
  
  
  name.for.loop <- paste0(virus.for.loop, "_cpo_weights")
  saving_path <-  paste0(
    "hhv_france/clean_data/output/sensitivity_analyses/weights/",
    name.for.loop,
    ".RDS"
  )
  if (file.exists(saving_path) == FALSE) {
    print("Changing to cpo")
    x <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      type_log_score = "loocv",
      input_list_weights = bma_draw$list_weights
    )
    
    saveRDS(x, saving_path)
    rm("x")
  } else{
    print("Already available")
  }
  
  
  name.for.loop <- paste0(virus.for.loop, "_gcpo25_weights")
  saving_path <-  paste0(
    "hhv_france/clean_data/output/sensitivity_analyses/weights/",
    name.for.loop,
    ".RDS"
  )
  if (file.exists(saving_path) == FALSE) {
    print("Changing to gcpo")
    x <- summary_analyses(
      virus.loop = virus.for.loop,
      n_draws = n_draws_choice,
      type_log_score = "lgocv.m25",
      input_list_weights = bma_draw$list_weights
    )
    saveRDS(x, saving_path)
    rm("x")
  } else{
    print("Already available")
  }
  
  assign(paste0(virus.for.loop, "_bma_analysis"),
         bma_draw)
  rm("bma_draw")
  gc(full = T)
}




# Saving image ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gc()
rm(list = setdiff(
  ls(),
  c(
    paste0(c("hsv1", "hsv2", "cmv", "vzv", "ebv"), "_bma_analysis"),
    "list_strata",
    "summary_analyses",
    "c_one_of"
  )
))
save.image("hhv_france/clean_data/output/results_bma.rda")
