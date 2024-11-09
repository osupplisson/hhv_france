# This code allows to fit all the models for all viruses
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")

# Import data ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
list2env(data_baseline, envir = .GlobalEnv)

# INLA OPTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting the maximum number of thread that INLA is able to use -----------------------------------------------------------------------------------------------------------------------------------
## set number of threads to 1 for ensure reproducibility
inla.setOption("num.threads", "12:2")

# For loop for fit --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set_virus <- c("cmv", "ebv", "vzv", "hsv1", "hsv2")


var_to_prec <- function(sigma) {
  tibble("var" = sigma ^ 2,
         "prec" = 1 / sigma ^ 2)
}
prec_to_var <- function(prec) {
  tibble("var" = 1 / prec,
         "std" = sqrt(1 / prec))
}

# ggplot(data.frame(x = fx(seq(-10,10), 5)), aes(x)) +
#   geom_density(colour = "blue")


# Loading list of model if it exists
try(set_formula <-
      readRDS("hhv_france/clean_data/output/status_model.RDS"))


if (!exists("set_formula")) {
  # If previous try failed: creating the tibble with formula
  set_formula <- expand.grid(
    age = c("aar1", "arw2"),
    spatial = c(
      "nb5",
      "nb10",
      "nb15",
      "nb20",
      "soi",
      "queen",
      "delaunay",
      "gabriel",
      "relative"
    ),
    model = c("baseline"),
    virus_name = set_virus
  ) %>%
    # Formula
    mutate(
      formula = paste(spatial, age, sep = "_"),
      model_name = paste("fit", virus_name, formula, model, sep = "_")
    ) %>%
    # Running status
    mutate(status = NA) %>% 
    arrange(desc(age)) 
}


vector_overseas <- c(95, 96, 97, 98, 99)


## https://link.springer.com/article/10.1007/s10654-008-9230-x

for (virus in set_virus) {
  data.for.fit <- data[[paste0("results_", virus)]] %>%
    mutate(
      sex_dummy = sex_numeric - 1,
      Intercept = 1,
      group_sex = ifelse(sex == "Male", 1, 2),
      group_area = ifelse(district_numeric %in% vector_overseas,
                          1,
                          2),
      district_numeric_female = ifelse(sex == "Male", NA, district_numeric),
      age_numeric_metropolitan = ifelse(!district_numeric %in% vector_overseas, age_numeric, NA),
      age_numeric_overseas = ifelse(district_numeric %in% vector_overseas, age_numeric, NA),
      district_numeric_metropolitan = ifelse(!district_numeric %in% vector_overseas, district_numeric, NA),
      district_numeric_metropolitan = ifelse(
        district_numeric_metropolitan == 100,
        95,
        district_numeric_metropolitan
      ),
      district_numeric_overseas = ifelse(district_numeric %in% vector_overseas, district_numeric, NA),
      district_numeric_overseas = district_numeric_overseas - 94,
      district_numeric_bis = district_numeric
    )
  
  
  
  set_formula_virus <- set_formula %>%
    filter(virus_name == virus)
  
  # Replace the set if it exists
  test <- NULL
  try(test <-
        readRDS(paste0(
          "hhv_france/clean_data/output/status_model_",
          virus,
          ".RDS"
        )))
  print(test)
  if (!is.null(test)) {
    set_formula_virus <- test
    rm("test")
  }
  
  for (row in seq(1, nrow(set_formula_virus %>% filter(virus_name == virus)))) {
    print(paste0(row, "/", nrow(
      set_formula_virus %>% filter(virus_name == virus)
    )))
    # Formula
    type_formula <- set_formula_virus %>%
      filter(virus_name == virus) %>%
      filter(row_number() == row) %>%
      .$formula
    
    
    # Model name
    name_model <- set_formula_virus %>%
      filter(virus_name == virus) %>%
      filter(row_number() == row) %>%
      .$model_name
    
    path_for_saving_fit <- paste0(path_to_fit, name_model)
    print(type_formula)
    print(path_for_saving_fit)
    
    statut_tmp <- set_formula_virus %>%
      filter(virus_name == virus) %>%
      filter(row_number() == row) %>%
      .$status
    
    if (is.na(statut_tmp)) {
      if (file.exists(paste0(path_for_saving_fit, ".RDS"))) {
        print("Object already available, import to check status")
        test <- readRDS(paste0(path_for_saving_fit, ".RDS"))
        print(test)
        if (class(test) == "inla") {
          statut_tmp <- "success"
        } else {
          statut_tmp <- "failure"
        }
        print(statut_tmp)
        set_formula_virus <- set_formula_virus %>%
          mutate(status = ifelse(
            virus_name == virus &
              row_number() == row,
            statut_tmp,
            status
          ))
        rm("test")
        rm("statut_tmp")
      }
    }
    
    # Model status
    status_model <- set_formula_virus %>%
      filter(virus_name == virus) %>%
      filter(row_number() == row) %>%
      .$status
    
    print(status_model)
    # If status is NA then no attemp to fit the model has been made so far
    if (is.na(status_model)) {
      print(virus)
      print(type_formula)
      print(name_model)
      print(path_to_fit)
      print(path_for_saving_fit)
      
      
      type_graph <-
        str_split(type_formula, pattern = "\\_")[[1]][1]
      if (str_detect(type_graph, "queen")) {
        graph_bym2_district <- inla_graph_districts_queen
      } else if (str_detect(type_graph, "delaunay")) {
        graph_bym2_district <- inla_graph_districts_delaunay
      } else if (str_detect(type_graph, "soi")) {
        graph_bym2_district <- inla_graph_districts_soi
      } else if (str_detect(type_graph, "gabriel")) {
        graph_bym2_district <- inla_graph_districts_gabriel
      } else if (str_detect(type_graph, "relative")) {
        graph_bym2_district <- inla_graph_districts_relative
      } else if (str_detect(type_graph, "nb")) {
        graph_bym2_district <-
          eval(parse(text = paste0(
            "inla_graph_districts_",
            type_graph
          )))
      }
      
      
      fx <- function(x) {
        2 * exp(x) / (1 + exp(x)) - 1
      }
      
      ggplot(data.frame(x = fx(rnorm(
        10000, 0, prec_to_var(0.1)$std
      ))), aes(x)) +
        geom_density(colour = "blue")
      
      
      district_effect <- paste0(
        'f(district_numeric,
                group = year_numeric,
                ngroup = 5,
        control.group = list(model = "ar1",
                hyper = list(rho = list(
                prior = "pc.cor1",
                param = c(0.5, 0.75)
                ))),',
        if (type_graph == "iid") {
          'model = "iid",
          hyper = list(
                  prec = list(
                    prior = "pc.prec",
                    param = c(u = 1,a = 0.01)
                  )
                ))'
        } else{
          'model = "bym2",
                graph = graph_bym2_district,
                adjust.for.con.comp = TRUE,
                scale.model = TRUE,
                hyper = list(
                  prec = list(
                    prior = "pc.prec",
                    param = c(u = 1,a = 0.01)
                  ),
                  phi = list(
                    prior = "pc",
                    param = c(
                      u = 0.5,
                      a = 0.5
                    )
                  )
                ))'
        }
      )
      
      
      # Choice for age
      if (str_detect(type_formula, "arw1")) {
        latent_age <- "rw1"
      } else if (str_detect(type_formula, "arw2")) {
        latent_age <- "rw2"
      } else if (str_detect(type_formula, "aar1")) {
        latent_age <- "ar1"
      } else if (str_detect(type_formula, "abym2")) {
        latent_age <- "bym2"
      } else if (str_detect(type_formula, "aiid")) {
        latent_age <- "iid"
      }
      
      
      fx <- function(x, n) {
        (exp(x) - 1) / (exp(x) + n - 1)
      }
      
      ggplot(data.frame(x = fx(
        rnorm(10000, 0, prec_to_var(0.2)$std), 2
      )), aes(x)) +
        geom_density(colour = "blue")
      
      
      
      age_effect <- paste0(
        'f(age_numeric,
        group = group_sex,
        ngroup = 2,
        control.group = list(model = "exchangeable",
        hyper = list(rho = list(prior = "normal",param = c(0,0.2)))),
      model = "',
      latent_age,
      '",',
      if (latent_age %in% c("rw1", "rw2", "bym2")) {
        "scale.model = TRUE,"
      },
      # "replicate = id.age.rep,",
      if (latent_age %in% c("rw1", "rw2", "iid")) {
        "hyper = list(
                       prec = list(
                         prior = 'pc.prec',
                         param = c(1, 0.01)
                       )
                     )"
      } else if (latent_age == "ar1") {
        "hyper = list(
                    prec = list(
                        prior = 'pc.prec',
                        param = c(1, 0.01)
                      ),
                      rho = list(
                        prior = 'pc.cor1',
                        param = c(0.5, 0.75)
                      )
                    )"
      } else if (latent_age == "bym2") {
        "graph = agraph,
                  hyper = list(
                    prec = list(
                        prior = 'pc.prec',
                        param = c(1, 0.01)
                      ),
                      phi = list(
                        prior = 'pc',
                        param = c(0.5, 0.5)
                      )
                    )"
      },
      ")"
      )
      
      
      
      
      # ## Females specific effect
      district_effect_female <-
        str_replace(district_effect,
                    "district_numeric",
                    "district_numeric_female")
      
      age_effect_metropolitan <- str_replace(
        age_effect,
        "age_numeric,",
        "age_numeric_metropolitan,
        replicate = district_numeric_metropolitan,"
      )
      
      age_effect_overseas <- str_replace(
        age_effect,
        "age_numeric,",
        "age_numeric_overseas,
        replicate = district_numeric_overseas,"
      )
      
      
      
      
      
      baseline_formula <- reformulate(
        c(
          "-1 + Intercept",
          district_effect,
          district_effect_female,
          age_effect,
          age_effect_metropolitan,
          age_effect_overseas
        ),
        "outcome"
      )
      
      
      if (latent_age == "rw2") {
        baseline_formula <- reformulate(
          c(
            "-1 + Intercept",
            district_effect,
            district_effect_female,
            age_effect
          ),
          "outcome"
        )
      }
      
      print(baseline_formula)
      # Prior for intercept
      list_prior_fixed <- list(
        mean = 0,
        mean.intercept = 0,
        prec = 1,
        prec.intercept = 1
      )
      
      # First run with gaussian/eb
      # Max running in second: 60*60=>1hour, most models should run in ~1/2 mins with this sample size
      inla.setOption(inla.timeout = 60 * 15)
      
      control_compute_inla <- list(
        openmp.strategy = "huge",
        dic = TRUE,
        waic = TRUE,
        config = TRUE
      )
      
      control_inla <- list(
        strategy = "laplace",
        int.strategy = "grid",
        optimise.strategy = "plain",
        fast = FALSE,
        stencil = 9,
        dz = 0.1,
        diff.logdens = 0.1,
        npoints = 100,
        numint.maxfeval = 80000000,
        stupid.search = T,
        optimiser = "gsl",
        tolerance = 0.0001
      )
      
      inla.pardiso.check()
      fit <- NULL
      fit <- try(inla(
        formula = baseline_formula,
        data = data.for.fit,
        family = "binomial",
        Ntrials = N,
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = control_compute_inla,
        control.inla = control_inla,
        control.fixed = list_prior_fixed,
        verbose = F,
        safe = T,
        debug = F
      ))
      
      
      
      print(summary(fit))
      # If fatal error during the initial fit: pass
      # Else: assigning NULL to the object, leading to error
      if (class(fit) == "try-error" | is.null(fit)) {
        print("Fatal error during the fit, pass")
        fit <- NULL
      } else {
        print(fit$misc$warnings)
      }
      
      if (!is.null(fit)) {
        if (function_check(fit) == T) {
          # If not null:
          ## Trying to properly refit everything up until: it is OK, or > max_rerun fit or fatal error
          rerun <- 1
          to_run <- T
          max_run <- 6
          while (to_run == T & rerun <= max_run) {
            if (rerun > 3) {
              print("Trying with GSL optimiser")
              fit$.args$control.inla$optimiser <- "gsl"
              fit$.args$control.inla$tolerance <- 0.0000001
            }
            print("Probably some issues in the fit, rerun")
            fit <- try(inla.rerun(fit))
            print(fit)
            if (class(fit) != "try-error") {
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
              fit <- NULL
              print("FAILURE :(")
            }
            
            if (rerun > max_run) {
              print("Number of max rerun reached")
              to_run <- F
              fit <- NULL
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
      
      if (!is.null(fit)) {
        fit$data_fit <- data.for.fit
        fit$fomula_good_format <- baseline_formula
        print("Computing CV")
        fit$loocv <-
          inla.group.cv(fit, num.level.sets = -1, strategy = "posterior")
        for (x in c(15, 25)) {
          print(x)
          namex <- paste0("lgocv.m", x)
          fit[[namex]] <-
            inla.group.cv(
              fit,
              num.level.sets = x,
              strategy = "posterior",
              size.max = x
            )
        }
      }
      
      
      # Function used for updating the status
      update_status <- function(df = set_formula_virus,
                                formula_input = type_formula,
                                virus_input = virus,
                                update = "success") {
        # Adding success
        df %>%
          mutate(status = ifelse(
            formula == formula_input &
              virus_name == virus_input,
            update,
            status
          ))
      }
      
      
      if (is.null(fit)) {
        set_formula_virus <- update_status(update = "failure")
      } else {
        set_formula_virus <- update_status(update = "success")
      }
      
      list_fit <- base::list.files(
        path_to_fit,
        pattern = ".RDS",
        all.files = T,
        full.names = FALSE
      )
      
      # Update status
      print(set_formula_virus %>%
              filter(virus_name == virus &
                       formula == type_formula))
      
      saveRDS(fit,
              paste0(path_for_saving_fit, ".RDS"))
      
      # Make sure it's deleted
      rm("fit")
    } else {
      print(name_model)
      print("Already estimated")
    }
  }
  
  # Updating the saving
  saveRDS(
    set_formula_virus,
    paste0(
      "hhv_france/clean_data/output/status_model_",
      virus,
      ".RDS"
    )
  )
}
