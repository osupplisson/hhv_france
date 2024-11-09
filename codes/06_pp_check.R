#This code performs the posterior predictive check
# Importing packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")


# Load data ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)

load("hhv_france/clean_data/output/results_bma.rda")
spatial_objects$list_dpt$geometry <- NULL

# Extracting hyperparameters values ---------------------------------------------------------------------------------------------------------------------------------------------------------------
#@ results_bma: One of the list produced by 04_bma.R
pp_check_functions_draws <- function(results_bma) {
  print("Importing models")
  
  # Vector to list
  list_model <- purrr::map(results_bma$list_model_to_import$model,
                           function(x) {
                             x
                           }) %>%
    # Set names
    purrr::set_names()
  
  
  # Import
  list_model <- list_model %>%
    purrr::map(function(x) {
      readRDS(str_replace_all(paste0(path_to_fit, x), "/", "//"))
    })
  
  
  print("PP checks")
  # Perform draws from the posterior joint for all models -----------------------------------------------------------------------------------------------------------------------------------------------------------
  # @set_fit: a list of inla fit
  # @n_sample: number of draws for each fit
  output <- purrr::map(names(list_model),
                       function(model.loop) {
                         print("Model to draw from:")
                         print(paste0(which(model.loop == names(list_model)), "/", length(list_model)))
                         print("Extract fit")
                         # Extract the fit
                         tmp <- list_model[[model.loop]]
                         # Extract the data used for the fit, keep only rows corresponding to the strata
                         data_used_for_fit <-
                           tmp$data_fit %>%
                           dplyr::mutate(id_row = dplyr::row_number()) %>%
                           dplyr::filter(name == "observed")
                         
                         print("Number of draws for the model")
                         # Extract the number of draws to perform for each model
                         n_to_draw <-
                           results_bma$weights %>%
                           dplyr::filter(model == model.loop) %>%
                           .$n_to_draw
                         
                         print("Sampling")
                         # Draws from the joint posterior
                         draw_approx_joint <-
                           inla.posterior.sample(n = n_to_draw,
                                                 result = tmp,
                                                 seed = 123)
                         
                         # Function to extract the parameter of the binomial
                         sim.p <- function() {
                           plogis(Predictor)
                         }
                         print("Extracting mu")
                         # Extract the parameter of the binomial
                         mu_draws <-
                           inla.posterior.sample.eval(sim.p,
                                                      draw_approx_joint)
                         print("Keeping only observed data")
                         # Keeping only rows from the strata in the fit
                         mu_draws <-
                           mu_draws[data_used_for_fit$id_row, ] %>%
                           as_tibble() %>%
                           tidytable()
                         
                         if (n_to_draw == 1) {
                           # If n_to_draw=1, changing the colname for consistency with n_to_draw>1
                           colnames(mu_draws) <- "sample:1"
                         }
                         print("Joining mu's with observed values")
                         # Joining
                         data_draws <- mu_draws %>%
                           cbind(data_used_for_fit) %>%
                           tidytable() %>%
                           tidytable::pivot_longer(
                             cols = starts_with("sample:"),
                             names_to = "sample",
                             values_to = "p"
                           ) %>%
                           tidytable::mutate(model = model.loop) %>%
                           virus_extraction(., "model") %>%
                           tidytable::mutate(all = "all") %>%
                           tidytable::select(
                             virus,
                             model,
                             all,
                             sex_numeric,
                             age_numeric,
                             district_numeric,
                             year_numeric,
                             N,
                             sample,
                             p,
                             outcome
                           ) %>%
                           left_join(
                             spatial_objects$list_dpt %>%
                               tidytable() %>%
                               tidytable::select(district_numeric, nom_region),
                             by = "district_numeric"
                           ) %>%
                           function_france_area()
                         print("Add draws")
                         ## Adding draws
                         data_draws <- data_draws %>%
                           tidytable::rowwise() %>%
                           tidytable::mutate(outcome_sim = rbinom(1,
                                                                  size = N,
                                                                  prob = p)) %>%
                           tidytable::ungroup()
                       }) %>%
    # Binding everything
    do.call(rbind, .) %>%
    #From tidytable to tibble
    tibble() %>%
    #Joining information
    left_join(age_vector, by = "age_numeric") %>%
    left_join(sex_vector, by = "sex_numeric") %>%
    left_join(year_vector, by = "year_numeric") %>%
    left_join(district_vector, by = "district_numeric") %>%
    function_france_area()
  
  return(output)
}



# Applying the function to all viruses --------------------------------------------------------------------------------------------------------------------------------------------------------------
output_saved <- base::list.files(
  "hhv_france/clean_data/output/post_fit_analyses/",
  pattern = ".RDS",
  all.files = T,
  full.names = FALSE
)

for (virus.for.loop in c("hsv1", "hsv2", "vzv", "ebv", "cmv")) {
  #Name for the pp
  name.for.loop <- paste0("pp_check_", virus.for.loop)
  #Name for the bma
  bma.for.loop <- paste0(virus.for.loop, "_bma_analysis")
  
  # If not available, perform PP check
  if (!paste0(name.for.loop, ".RDS") %in% output_saved) {
    print("Perform PP check for:")
    print(virus.for.loop)
    
    #Function
    #We are going to use tidytable here, to speed thing up
    #The bottleneck in terms of computation is the last part where we draw per rows
    library("tidytable")
    pp_check <-
      pp_check_functions_draws(results_bma = eval(parse(text = bma.for.loop)))
    #Summarising - Tidytable should be remove for group_at() function to work properly
    detach("package:tidytable", unload = TRUE)
    
    ##Producing figures
    graph_ppcheck <- function(df_input = NA ,
                              var_to_stratify = "virus",
                              var_to_stratify_2 = NULL,
                              type_plot = "ecdf",
                              labs_x_proportion = element_blank()) {
      removing_df <- function(plot) {
        plot <- ggplot2::ggplotGrob(plot)
        ggpubr::as_ggplot(plot)
      }
      
      df_input <- df_input %>%
        mutate(virus = str_to_upper(virus),
               virus = str_replace_all(virus, "HSV", "HSV-"))
      
      if (type_plot == "ecdf") {
        p <- df_input %>%
          ggplot()  +
          stat_ecdf(aes(x = 100 * outcome / N,
                        color = "A"), geom = "step") +
          stat_ecdf(aes(x = 100 * outcome_sim / N,
                        color = "B"), geom = "step") +
          labs(y = "Cumulative distribution",
               x = "Proportion of positive tests (in %)",
               color = element_blank()) +
          theme_classic() +
          scale_color_manual(
            values = c("darkred", "darkblue"),
            label = c("Observed", "Posterior")
          ) +
          theme(legend.position = "bottom",
                legend.direction = "horizontal") +
          guides(colour = guide_legend(title.position = "top",
                                       nrow = 1))
      } else if (type_plot == "density") {
        p <- df_input %>%
          ggplot()  +
          geom_histogram(
            aes(
              y = ..density..,
              x = 100 * outcome / N,
              fill = "A",
              color = "A"
            ),
            alpha = 0.1,
            binwidth = 2.5
          ) +
          geom_histogram(
            aes(
              y = ..density..,
              x = 100 * outcome_sim / N,
              fill = "B",
              color = "B"
            ),
            alpha = 0.1,
            binwidth = 2.5
          )  +
          labs(
            y = "Density",
            x = "Proportion of positive tests (in %)",
            fill = element_blank(),
            color = element_blank()
          ) +
          theme_classic() +
          scale_fill_manual(
            values = c("darkred", "darkblue"),
            label = c("Observed", "Posterior")
          ) +
          scale_color_manual(
            values = c("darkred", "darkblue"),
            label = c("Observed", "Posterior")
          ) +
          theme(legend.position = "bottom",
                legend.direction = "horizontal") +
          guides(colour = guide_legend(title.position = "top",
                                       nrow = 1))
      } else if (type_plot == "proportion") {
        df_input <- rbind(
            df_input %>% mutate(y = 100 * outcome / N, name = "observed"),
            df_input %>% mutate(y = 100 * outcome_sim / N, name = "simulated")
          ) %>%
          mutate(bucket = base::cut(
            y,
            breaks = seq(0, 100, 5),
            include.lowest = T
          )) %>%
          group_by_at(
            vars(
              "sample",
              "model",
              "bucket",
              "name",
              var_to_stratify,
              var_to_stratify_2
            )
          ) %>%
          count() %>%
          ungroup() %>%
          group_by_at(vars(
            "sample",
            "model",
            "name",
            var_to_stratify,
            var_to_stratify_2
          )) %>%
          mutate(pct = 100 * n / sum(n)) %>%
          dplyr::select(-n) %>%
          pivot_wider(names_from = name,
                      values_from = "pct") %>%
          ungroup()
        
        
        
        df_input <- df_input %>%
          group_by_at(vars("bucket", var_to_stratify, var_to_stratify_2)) %>%
          rename(value = simulated) %>%
          function_for_summary_draws() %>%
          ungroup() %>%
          left_join(
            df_input %>%
              group_by(bucket) %>%
              rename(value = observed) %>%
              function_for_summary_draws() %>%
              ungroup() %>%
              dplyr::select(bucket,
                            observed = mean),
            by = c('bucket')
          )
        
        p <- df_input %>%
          ggplot() +
          geom_point(aes(
            x = bucket,
            y = observed,
            color = "A"
          ),
          alpha = 0.3) +
          geom_point(aes(
            x = bucket,
            y = mean,
            color = "B"
          ),
          alpha = 0.3) +
          geom_errorbar(aes(
            x = bucket,
            ymin = qi_lb,
            ymax = qi_ub,
            color = "C"
          ),
          alpha = 0.3) +
          scale_color_manual(
            values = c("blue",
                       "darkred",
                       "black"),
            label = c("Observed",
                      "Posterior average",
                      "Posterior ETI95%")
          ) +
          labs(y = "% of positive tests",
               x = labs_x_proportion,
               color = element_blank()) +
          theme_classic() +
          theme(
            legend.position = "bottom",
            legend.box = "vertical",
            legend.margin = margin()
          ) +
          theme(axis.text.x = element_text(
            angle = -90,
            vjust = 0,
            hjust = 0
          ))
      } 
      if (type_plot != "calibration") {
        if (!is.null(var_to_stratify_2)) {
          p <-
            p + facet_grid(df_input[[var_to_stratify_2]] ~ df_input[[var_to_stratify]], scales = "free")
        } else{
          p <- p + facet_grid( ~ df_input[[var_to_stratify]], scales = "free")
        }
      }
      p <- removing_df(p)
      return(p)
    }
    print(pp_check)
    pp_graph <- list(
      "pp_all" =  graph_ppcheck(df_input = pp_check,
                                type_plot = "ecdf"),
      "pp_sex" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "sex",
        type_plot = "ecdf"
      ),
      "pp_area" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "france_area",
        type_plot = "ecdf"
      ),
      "pp_region" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "nom_region",
        type_plot = "ecdf"
      ),
      "pp_age" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "age_factor",
        type_plot = "ecdf"
      ),
      "pp_year" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "year",
        type_plot = "ecdf"
      ),
      "pp_all_density" =  graph_ppcheck(df_input = pp_check,
                                        type_plot = "density"),
      "pp_sex_density" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "sex",
        type_plot = "density"
      ),
      "pp_area_density" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "france_area",
        type_plot = "density"
      ),
      "pp_region_density" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "nom_region",
        type_plot = "density"
      ),
      "pp_age_density" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "age_factor",
        type_plot = "density"
      ),
      "pp_year_density" = graph_ppcheck(
        df_input = pp_check,
        var_to_stratify_2 = "year",
        type_plot = "density"
      ),
      "pp_all_proportion" =  graph_ppcheck(df_input = pp_check,
                                           type_plot = "proportion"),
      "pp_sex_proportion" =  graph_ppcheck(df_input = pp_check,
                                           var_to_stratify_2 = "sex",
                                           type_plot = "proportion")
    )
    rm("pp_check")
    gc(full = T)
    
    # Save
    saveRDS(
      pp_graph,
      paste0(
        "hhv_france/clean_data/output/post_fit_analyses/",
        name.for.loop,
        ".RDS"
      )
    )
    
    # Assign
    assign(name.for.loop,
           pp_graph)
    rm("pp_graph")
  } else {
    print("PP check already available for:")
    print(virus.for.loop)
    # Otherwise download
    assign(name.for.loop,
           readRDS(
             paste0(
               "hhv_france/clean_data/output/post_fit_analyses/",
               name.for.loop,
               ".RDS"
             )
           ))
  }
}




gc()
rm(list = setdiff(ls(), paste0(
  "pp_check_", c("hsv1", "hsv2", "cmv", "vzv", "ebv")
)))
save.image("hhv_france/clean_data/output/pp_check.rda")
