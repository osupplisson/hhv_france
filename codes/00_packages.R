###This codes import all the package needed for conducting the analysis
###It also defines some useful functions
# Cleaning ent----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Cleaning ent")
rm(list = ls())
options(scipen=999)

# Cleaning packages -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Cleaning packages")
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

# Packages -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("tidybayes")
library("tidyverse")
library("tidylog")
library("spdep")

print("MAKE SURE TO USE INLA LAST TESTING VERSION !!!!")
print("You may use the following code to get it:")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library("INLA")

# Seed --------------------------------------------------------------------
set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")

# Choice start data -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
start_obs <- "2018-01-01"

# Convenience function ----------------------------------------------------
function_for_summary_draws <- function(df) {
  df %>%
    dplyr::summarise(
      mean = mean(value, na.rm = T),
      median = median(value, na.rm = T),
      # mode = Mode(value),
      # hdi_lb = hdi(value)[1],
      # hdi_ub = hdi(value)[2],
      # hdci_lb = hdci(value)[1],
      # hdci_ub = hdci(value)[2],
      qi_lb = tidybayes::qi(value, na.rm = T)[1],
      qi_ub = tidybayes::qi(value, na.rm = T)[2],
      qi_lb_09 = tidybayes::qi(value, .width = 0.9, na.rm = T)[1],
      qi_ub_09 = tidybayes::qi(value, .width = 0.9, na.rm = T)[2],
      qi_lb_08 = tidybayes::qi(value, .width = 0.8, na.rm = T)[1],
      qi_ub_08 = tidybayes::qi(value, .width = 0.8, na.rm = T)[2],
      qi_lb_07 = tidybayes::qi(value, .width = 0.7, na.rm = T)[1],
      qi_ub_07 = tidybayes::qi(value, .width = 0.7, na.rm = T)[2],
      qi_lb_06 = tidybayes::qi(value, .width = 0.6, na.rm = T)[1],
      qi_ub_06 = tidybayes::qi(value, .width = 0.6, na.rm = T)[2],
      qi_lb_05 = tidybayes::qi(value, .width = 0.5, na.rm = T)[1],
      qi_ub_05 = tidybayes::qi(value, .width = 0.5, na.rm = T)[2],
      min = min(value, na.rm = T),
      max = max(value, na.rm = T),
      higher_0 = mean(value > 0, na.rm = T),
      higher_1 = mean(value > 1, na.rm = T),
      higher_05 = mean(value > 0.5, na.rm = T),
      lower_0 = mean(value < 0, na.rm = T),
      lower_minus1 = mean(value < 1, na.rm = T),
      lower_minus05 = mean(value < 0.5, na.rm = T)
    )
}



# Path --------------------------------------------------------------------
#These are some paths that are going to be useful in the process
path_to_inla_directory <- "hhv_france/inla_run/"
path_to_fit <- "hhv_france/clean_data/output/fit/"
path_to_post_fit_analyses <-
  "hhv_france/clean_data/output/post_fit_analyses/"



# Function for extracting virus name --------------------------------------------------------------------------------------------------------------------------------------------------------------
#These is a function useful to extract the virus name from a string
virus_extraction <- function(df,
                             variable) {
  df %>%
    mutate(
      virus = case_when(
        str_detect(!!!rlang::syms(variable), "hsv1") ~ "hsv1",
        str_detect(!!!rlang::syms(variable), "hsv2") ~ "hsv2",
        str_detect(!!!rlang::syms(variable), "cmv") ~ "cmv",
        str_detect(!!!rlang::syms(variable), "vzv") ~ "vzv",
        str_detect(!!!rlang::syms(variable), "ebv") &
          !str_detect(!!!rlang::syms(variable), "_eba") &
          !str_detect(!!!rlang::syms(variable), "_vca") ~ "ebv",
        str_detect(!!!rlang::syms(variable), "_eba") ~ "ebv (ebna)",
        str_detect(!!!rlang::syms(variable), "_vca") ~ "ebv (vca)"
      )
    )
}



# Extracting model name ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Extracting model name
extraction_model_name <- function(df) {
  df %>%
    mutate(
      label_age =  case_when(
        str_detect(model, "aiid") ~ "IID",
        str_detect(model, "arw1") ~ "RW(1)",
        str_detect(model, "arw2") ~ "RW(2)",
        str_detect(model, "aar1") ~ "AR(1)",
        str_detect(model, "abym2") ~ "BYM2",
        TRUE ~ ""
      ),
      label_time = case_when(
        str_detect(model, "tiid") ~ "IID",
        str_detect(model, "trw1") ~ "RW(1)",
        str_detect(model, "trw2") ~ "RW(2)",
        str_detect(model, "tar1") ~ "AR(1)",
        str_detect(model, "tbym2") ~ "BYM2",
        TRUE ~ ""
      ),
      label_district = case_when(
        str_detect(model, "diid") ~ "IID",
        str_detect(model, "queen") ~ "BYM2(Queen)",
        str_detect(model, "delaunay") ~ "BYM2(Delaunay)",
        str_detect(model, "soi") ~ "BYM2(SOI)",
        str_detect(model, "gabriel") ~ "BYM2(Gabriel)",
        str_detect(model, "relative") ~ "BYM2(Relative)",
        str_detect(model, "nb5") ~ "BYM2(NN5)",
        str_detect(model, "nb10") ~ "BYM2(NN10)",
        str_detect(model, "nb15") ~ "BYM2(NN15)",
        str_detect(model, "nb20") ~ "BYM2(NN20)",
        str_detect(model, "nb25") ~ "BYM2(NN25)",
        TRUE ~ ""
      ),
      label_sex = case_when(
        str_detect(model, "siid") ~ "IID",
        str_detect(model, "sbym2") ~ "BYM2",
        TRUE ~ 'Fixed'
      ),
      model_name = paste0(
        "Age:",
        label_age,
        "-Time:",
        label_time,
        "-Space:",
        label_district
      )
    )
}

#Ordering the age classes
ordering_age <- function(df) {
  df %>%
    mutate(age_factor = factor(
      age_class,
      levels = c(
        "[0,4]",
        "(4,9]",
        "(9,14]",
        "(14,19]",
        "(19,24]",
        "(24,29]",
        "(29,34]",
        "(34,39]",
        "(39,44]",
        "(44,49]",
        "(49,54]",
        "(54,59]",
        "(59,64]",
        "(64,69]",
        "(69,74]",
        "(74,79]",
        "(79,84]",
        "(84,89]",
        "(89,94]",
        "95+"
      )
    ))
}

# Function for ordering virus
order_virus <- function(df) {
  df %>%
    mutate(
      virus = case_when(
        virus == "hsv1" ~ "HSV-1",
        virus == "hsv2" ~ "HSV-2",
        virus == "cmv" ~ "CMV",
        virus == "vzv" ~ "VZV",
        virus == "ebv" ~ "EBV",
        virus == "ebv (vca)" ~ "EBV (VCA)",
        virus == "ebv (ebna)" ~ "EBV (EBNA)"
      ),
      virus = factor(
        virus,
        levels = c("HSV-1",
                   "HSV-2",
                   "VZV",
                   "EBV",
                   "EBV (EBNA)",
                   "EBV (VCA)",
                   "CMV")
      )
    )
}


# France are
#Function to create a new variables with huge areas, starting from region
function_france_area <- function(df) {
  df %>%
    mutate(
      france_area = case_when(
        nom_region %in% c("Guadeloupe", "Martinique", "La Réunion", "Guyane", "Mayotte") ~ "Overseas France",
        TRUE ~ "Metropolitan France"
      )
    )
}



# Checking INLA -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function_check <- function(fit) {
  ## If refit worked: either it's ok and stop or we try to refit if rerun<10
  a <- fit$mode$mode.status
  b <- identical(fit$misc$warnings, character(0))
  c <- fit$waic$waic
  d <- fit$dic$dic
  e <-
    sum(abs(fit$misc$cov.intern[upper.tri(fit$misc$cov.intern)]))
  if (nrow(fit$misc$cov.intern) == 0) {
    # Because there is no cov.intern matrix in this case
    e <- 1
  }
  f <- fit$ok
  print("Mode status, should be 0:")
  print(a)
  print("Warning during fit, should be TRUE:")
  print(b)
  print("DIC and WAIC should be non infinite")
  print(c)
  print(d)
  print("Should be != 0:")
  print(e)
  print("Should be TRUE")
  print(f)
  print("Checking that everything is fine")
  if (a > 0 |
      b != T |
      is.infinite(c) |
      is.infinite(d) |
      e == 0 |
      f != TRUE) {
    # Issue: rerun
    T
  } else {
    # No issue
    F
  }
}






