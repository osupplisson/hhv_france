## This code builds the aggregated (shared) dataset
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")


# Defining functions ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function_data_for_fit <- function(type_analysis = "baseline") {
  # Import data ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (type_analysis == "baseline") {
    data <- readRDS("hhv_france/clean_data/input_models/data_clean.RDS")
  }else if(type_analysis == "betabinomial"){
    data <- readRDS("hhv_france/clean_data/input_models/data_clean_for_betabinomial.RDS")
  }
  
  #Spatial objects -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  spatial_objects <-
    readRDS("hhv_france/clean_data/input_models/spatial_objects.RDS")
  spatial_objects$list_dpt$geometry <- NULL
  
  # Strata ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  strata <- readRDS("hhv_france/clean_data/input_models/strata.RDS")
  
  
  
  # Year --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  data <- data %>%
    mutate(year = format(date_sampling, "%Y"))
  
  
  # Age class ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  age_cut <- c(unique(strata$lb) - 1, 150)
  age_cut[age_cut == -1] <- 0
  
  data <- data %>%
    mutate(age_class = cut(age, age_cut, include.lowest = TRUE))
  
  data <- data %>%
    mutate(
      age_class = as.character(age_class),
      age_class = ifelse(age_class == "(94,150]", "95+", age_class)
    )
  
  # Sex ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  data <- data %>%
    mutate(sex_numeric = ifelse(sex == "Male",
                                1,
                                2))
  table(data$sex_numeric)
  # Adding id ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  unique_year <- strata %>% distinct(year, year_numeric, year_factor)
  unique_age <-
    strata %>% distinct(age_class, age_numeric, age_factor)
  
  data <- data %>%
    left_join(unique_year, by = "year") %>%
    left_join(unique_age, by = "age_class")
  
  
  # Adding region names -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  tmp <- spatial_objects$list_dpt
  tmp$geometry <- NULL
  data <- data %>%
    left_join(tmp %>%  select(nom, nom_region), by = 'nom')
  rm('tmp')
  
  # INLA GRAPHS -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  # IMPORTING ALL INLA Graphs -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  path_to_graph_inla <- "hhv_france/clean_data/input_models/"
  
  list_inla_graph <- base::list.files(
    path_to_graph_inla,
    pattern = ".graph",
    all.files = T,
    full.names = FALSE
  )
  
  for (graph.loop in list_inla_graph) {
    print(graph.loop)
    name_object <- str_remove(graph.loop, "\\.graph")
    link_object <- paste0(path_to_graph_inla, graph.loop)
    graph <- inla.read.graph(link_object)
    assign(name_object,
           graph)
  }
  
  
  # Creating graph for time -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  tgraph <-
    sparseMatrix(i = c(2:max(unique(
      data$year_numeric
    )), 1:(max(
      unique(data$year_numeric)
    ) - 1)),
    j = c(1:(max(
      unique(data$year_numeric)
    ) - 1), 2:max(unique(
      data$year_numeric
    ))),
    x = 1)
  
  
  # Creating graph for age -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  agraph <- sparseMatrix(i = c(2:max(unique(data$age_numeric)), 1:(max(
    unique(data$age_numeric)
  ) - 1)),
  j = c(1:(max(
    unique(data$age_numeric)
  ) - 1), 2:max(unique(data$age_numeric))),
  x = 1)
  
  # Several vectors that may be useful --------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  sex_vector <- strata %>%
    dplyr::select(sex_numeric,
                  sex = sex,
                  sex_factor = sex_factor) %>%
    distinct(sex, .keep_all = T)
  
  
  age_vector <- strata %>%
    dplyr::select(age_numeric,
                  age_class) %>%
    distinct(age_class, .keep_all = T) %>%
    ordering_age()
  
  
  district_vector <- strata %>%
    dplyr::select(district_numeric,
                  nom) %>%
    distinct(nom, .keep_all = T)
  
  
  year_vector <- strata %>%
    dplyr::select(year,
                  year_numeric) %>%
    distinct(year, .keep_all = T)
  
  
  
  
  # Cleaning ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  data <- data %>%
    select(
      id_row_init,
      id_patient,
      sex,
      sex_numeric,
      age_class,
      age_numeric,
      age_factor,
      nom,
      id_district = district,
      district_numeric,
      centroid_long,
      centroid_lat,
      year,
      year_numeric,
      year_factor,
      name,
      value
    )
  
  # saving image ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  rm("function_for_summary_draws")
  rm("path_to_inla_directory")
  rm("path_to_fit")
  rm("path_to_post_fit_analyses")
  rm("extraction_model_name")
  rm("ordering_age")
  rm("order_virus")
  rm("function_france_area")
  
  # Creating unique id by subgroups --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  replicate <- expand.grid(
    sex_numeric = unique(sex_vector$sex_numeric),
    age_numeric = unique(age_vector$age_numeric),
    year_numeric = unique(year_vector$year_numeric),
    district_numeric = unique(district_vector$district_numeric)
  ) %>%
    group_by(age_numeric, year_numeric, district_numeric) %>%
    mutate(id.sex.rep = row_number()) %>%
    group_by(sex_numeric, year_numeric, district_numeric) %>%
    mutate(id.age.rep = row_number()) %>%
    group_by(sex_numeric, age_numeric, district_numeric) %>%
    mutate(id.year.rep = row_number()) %>%
    group_by(sex_numeric, age_numeric, year_numeric) %>%
    mutate(id.district.rep = row_number()) %>%
    group_by(sex_numeric, age_numeric, year_numeric, district_numeric) %>%
    mutate(id.strata = row_number()) %>%
    ungroup()
  
  
  # Data used for fit ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  keys <- data %>% group_by(name) %>% group_keys()
  data_for_fit <- data %>%
    virus_extraction("name") %>%
    group_split(name, .keep = TRUE) %>%
    purrr::imap(
      ~ .x %>%
        #Numeric outcome
        mutate(outcome = ifelse(value == "p", 1, 0)) %>%
        #Aggregation by stratum
        group_by(
          sex_numeric,
          age_numeric,
          year_numeric,
          district_numeric,
          virus
        ) %>%
        summarise(N = n(),
                  outcome = sum(outcome)) %>%
        ungroup()
    )
  
  
  #Checking that everything went fine
  data %>%
    group_by(name) %>%
    count()
  
  sum(data_for_fit[[1]]$N)
  sum(data_for_fit[[2]]$N)
  sum(data_for_fit[[3]]$N)
  sum(data_for_fit[[4]]$N)
  
  data_for_fit_final <- data_for_fit %>%
    purrr::imap(
      #Adding all the strata for prediction to the actual grid
      ~ .x %>%
        mutate(name = "observed") %>%
        plyr::rbind.fill(
          strata %>%
            mutate(name = "predicted",
                   outcome = NA) %>%
            dplyr::select(one_of(colnames(.x), "name"))
        ) %>%
        #Adding id
        mutate(id_data = row_number()) %>%
        #Adding replicate id if needed
        left_join(
          replicate,
          by = c(
            "sex_numeric",
            "age_numeric",
            "year_numeric",
            "district_numeric"
          )
        ) %>%
        #Adding label
        left_join(
          strata %>% distinct(sex_numeric, .keep_all = T) %>% dplyr::select(sex, sex_numeric, sex_factor),
          by = 'sex_numeric'
        ) %>%
        left_join(
          strata %>% distinct(year_numeric, .keep_all = T) %>% dplyr::select(year, year_numeric, year_factor),
          by = 'year_numeric'
        ) %>%
        left_join(
          strata %>% distinct(age_numeric, .keep_all = T) %>% dplyr::select(age_class, age_numeric, age_factor),
          by = 'age_numeric'
        ) %>%
        left_join(
          strata %>% distinct(district_numeric, .keep_all = T) %>% dplyr::select(
            nom,
            centroid_long,
            centroid_lat,
            id_district,
            district_numeric
          ),
          by = 'district_numeric'
        ) %>%
        left_join(
          spatial_objects$list_dpt  %>% dplyr::select(nom_region, district_numeric),
          by = 'district_numeric'
        )
    )
  
  
  
  sum(data_for_fit[[1]]$N)
  sum(data_for_fit_final[[1]] %>% filter(name == "observed") %>% .$N)
  sum(data_for_fit[[2]]$N)
  sum(data_for_fit_final[[2]] %>% filter(name == "observed") %>% .$N)
  sum(data_for_fit[[3]]$N)
  sum(data_for_fit_final[[3]] %>% filter(name == "observed") %>% .$N)
  sum(data_for_fit[[4]]$N)
  sum(data_for_fit_final[[4]] %>% filter(name == "observed") %>% .$N)
  
  
  #Adding keys
  names(data_for_fit_final) <- keys$name
  
  
  table(data_for_fit_final$results_cmv %>% filter(name == "observed") %>% .$year)
  table(data_for_fit_final$results_ebv %>% filter(name == "observed") %>% .$year)
  table(data_for_fit_final$results_vzv %>% filter(name == "observed") %>% .$year)
  table(data_for_fit_final$results_hsv1 %>% filter(name == "observed") %>% .$year)
  table(data_for_fit_final$results_hsv2 %>% filter(name == "observed") %>% .$year)
  table(data_for_fit_final$results_hsv2 %>% filter(name == "predicted") %>% .$year)
  
  data <- data_for_fit_final
  
  rm("data_for_fit")
  rm("data_for_fit_final")
  rm("keys")
  rm("virus_extraction")
  
  # Reloading spatial object to get the geometry right ----------------------------------------------------------------------------------------------------------------------------------------------
  spatial_objects <-
    readRDS("hhv_france/clean_data/input_models/spatial_objects.RDS")
  
  return(as.list(environment()))
}

data_baseline <- function_data_for_fit(type_analysis = "baseline")
data_betabinomial <- function_data_for_fit(type_analysis = "betabinomial")
save.image("hhv_france/clean_data/input_models/all_dataset_ready_for_fit.rda")
