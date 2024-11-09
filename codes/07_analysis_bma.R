#This code builds a .RDA used in the .RMD and defines some functions used in the RMD
# Importing packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")


# Importing data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/output/results_bma.rda")
load("hhv_france/clean_data/output/pp_check.rda")
load("hhv_france/clean_data/output/post_fit_analyses/hyperparameters.rda")
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)




# Saving all model status -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

saveRDS(
  rbind(
    readRDS("hhv_france/clean_data/output/status_model_hsv1.RDS"),
    readRDS("hhv_france/clean_data/output/status_model_hsv2.RDS"),
    readRDS("hhv_france/clean_data/output/status_model_vzv.RDS"),
    readRDS("hhv_france/clean_data/output/status_model_cmv.RDS"),
    readRDS("hhv_france/clean_data/output/status_model_ebv.RDS")
  ),
  "hhv_france/clean_data/output/status_model.RDS"
)

# Binding the post-stratified results -------------------------------------------------------------------------------------------------------------------------------------------------------------
binding_strata <- function(type,
                           var) {
  rbind(
    hsv1_bma_analysis[[type]][[var]] %>%
      mutate(virus = "hsv1"),
    hsv2_bma_analysis[[type]][[var]] %>%
      mutate(virus = "hsv2"),
    cmv_bma_analysis[[type]][[var]] %>%
      mutate(virus = "cmv"),
    vzv_bma_analysis[[type]][[var]] %>%
      mutate(virus = "vzv"),
    ebv_bma_analysis[[type]][[var]] %>%
      mutate(virus = "ebv")
  ) %>%
    order_virus()
}

# Table weights -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
all_weights <- rbind(
  hsv1_bma_analysis$weights,
  hsv2_bma_analysis$weights,
  cmv_bma_analysis$weights,
  vzv_bma_analysis$weights,
  ebv_bma_analysis$weights
)



# # Table hyperparameters -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hyperparam <- rbind(
  hyperparameters_hsv1,
  hyperparameters_hsv2,
  hyperparameters_cmv,
  hyperparameters_vzv,
  hyperparameters_ebv
)


plot_hyperparam <-
  function(param_choice = "Phi for district_numeric",
           label_choice = NULL) {
    if (is.null(label_choice)) {
      label_choice <- str_remove(param_choice, "\\_numeric")
    }
    hyperparam %>%
      filter(param == param_choice) %>%
      order_virus() %>%
      ggplot() +
      geom_point(aes(
        y = fct_reorder(model_name, mean),
        x = mean,
        color = "A"
      )) +
      geom_errorbar(aes(
        y = fct_reorder(model_name, mean),
        xmin = `0.025quant`,
        xmax = `0.975quant`,
        color = "B"
      )) +
      labs(
        y = "Models",
        x = element_blank(),
        color = element_blank(),
        subtitle = label_choice
      ) +
      scale_color_manual(
        values = c("darkblue",
                   "darkred"),
        labels = c("Posterior average",
                   "Posterior ETI95%")
      ) +
      facet_wrap( ~ virus) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))
  }


# Table seroprev -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------<
faceted_plot_seroprev <- function(df_plot,
                                  var_to_stratify) {
  df_plot %>%
    ggplot() +
    geom_point(aes(x = mean,
                   y = .data[[var_to_stratify]],
                   color = "A")) +
    geom_errorbar(aes(
      width = 0.5,
      xmin = qi_lb,
      xmax = qi_ub,
      y = .data[[var_to_stratify]],
      color = "B",
      width = 0.5
    )) +
    geom_point(aes(x = observed,
                   y = .data[[var_to_stratify]],
                   color = "C")) +
    labs(y = element_blank(),
         x = element_blank(),
         color = element_blank()) +
    facet_wrap( ~ virus, nrow = 1, scale = "free_x") +
    scale_color_manual(
      values = c("blue", "black", "red"),
      labels = c("Posterior average",
                 "Posterior ETI95%",
                 "Observed")
    ) +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    scale_x_continuous(labels =  scales::comma,
                       n.breaks = 10) +
    theme(panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.1
    )) +
    scale_y_discrete(expand = waiver()) +
    theme(strip.placement = 'outside')
}
graph_seroprev <- function(var_input = "age_numeric",
                           var_to_stratify = "age_factor") {
  df_seroprev <- binding_strata(type = "bma_results_summary",
                                var = var_input)
  
  
  #Binding the data
  observed_data <- data  %>%
    do.call(rbind, .)  %>%
    order_virus() %>%
    filter(name == "observed")
  
  
  
  if (var_input %in% c("sex_numeric",
                       "c(\"france_area\", \"year_numeric\")")) {
    obs <- observed_data %>%
      rename(sex = sex) %>%
      dplyr::group_by(sex, virus) %>%
      dplyr::summarise(observed = 100 * sum(outcome) / sum(N))
  }
  
  if (var_input %in% c("age_numeric",
                       "france_area",
                       "nom_region",
                       "district_numeric",
                       "year_numeric")) {
    obs <- observed_data %>%
      function_france_area() %>%
      dplyr::group_by_at(c(var_input, "virus")) %>%
      dplyr::summarise(observed = 100 * sum(outcome) / sum(N))
  }
  
  if (str_detect(var_input, "diff")) {
    var_input <- case_when(
      str_detect(var_input, "sex") ~ "sex_numeric",
      str_detect(var_input, "age") ~ "age_numeric",
      str_detect(var_input, "year") ~ "year_numeric"
    )
    
    obs <- observed_data %>%
      dplyr::group_by_at(c(var_input, "virus"))  %>%
      dplyr::summarise(observed = 100 * sum(outcome) / sum(N)) %>%
      dplyr::group_by_at("virus") %>%
      arrange_at(var_input) %>%
      mutate(value = observed - lag(observed)) %>%
      dplyr::select(-observed) %>%
      rename(observed = value) %>%
      ungroup()
  } else if (str_detect(var_input, "special")) {
    var_input <- case_when(
      str_detect(var_input, "district") ~ "nom",
      str_detect(var_input, "region") ~ "nom_region"
    )
    
    
    obs <- observed_data %>%
      dplyr::group_by_at(c("virus", var_input)) %>%
      dplyr::summarise(observed = 100 * sum(outcome) / sum(N)) %>%
      ungroup()
    
    if (var_input  == "nom") {
      ref_df <- obs[obs$nom == "Paris",]
      ref_df <-
        ref_df %>% rename(observed_ref = observed, nom_ref = nom)
    } else {
      ref_df <- obs[obs$nom_region == "Île-de-France",]
      ref_df <-
        ref_df %>% rename(observed_ref = observed, nom_ref = nom_region)
    }
    
    obs <- obs %>%
      left_join(ref_df, by = 'virus') %>%
      mutate(value = observed - observed_ref) %>%
      dplyr::select(-observed) %>%
      rename(observed = value) %>%
      ungroup()
  }
  
  if (exists("obs")) {
    df_plot <- df_seroprev %>%
      left_join(obs) %>%
      filter(!is.na(mean))
  } else{
    df_plot <- df_seroprev %>%
      mutate(observed = NA)
  }
  
  if (var_input == "france_area" &
      var_to_stratify == "france_area") {
    df_seroprev <- binding_strata(type = "bma_results_summary",
                                  var = var_input)
    df_seroprev_all <- binding_strata(type = "bma_results_summary",
                                      var = "all")
    df_seroprev <- plyr::rbind.fill(df_seroprev,
                                    df_seroprev_all)  %>%
      mutate(
        france_area = case_when(
          str_detect(france_area, "Metropolitan") ~ "Metropolitan",
          str_detect(france_area, "Overseas") ~ "Overseas",
          TRUE ~ "France"
        )
      )
    
    obs <- observed_data %>%
      function_france_area() %>%
      dplyr::group_by_at(c(var_input, "virus")) %>%
      dplyr::summarise(observed = 100 * sum(outcome) / sum(N)) %>%
      plyr::rbind.fill(
        observed_data %>%
          dplyr::group_by_at("virus") %>%
          dplyr::summarise(observed = 100 * sum(outcome) / sum(N))
      ) %>%
      mutate(
        france_area = case_when(
          str_detect(france_area, "Metropolitan") ~ "Metropolitan",
          str_detect(france_area, "Overseas") ~ "Overseas",
          TRUE ~ "France"
        )
      )
    
    
    df_plot <- left_join(df_seroprev,
                         obs, by = c("virus", "france_area")) %>%
      mutate(france_area = factor(
        france_area,
        levels = c("Overseas",
                   "Metropolitan",
                   "France")
      ))
  }
  
  
  faceted_plot_seroprev(df_plot = df_plot,
                        var_to_stratify = var_to_stratify)
  
  
}

graph_seroprev(var_input = "sex_numeric",
               var_to_stratify = "sex")

save.image("hhv_france/clean_data/output/all_data_for_rmd.RDS")
