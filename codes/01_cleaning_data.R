#This code: importing and cleaning the data
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")

# Importing data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <-
  readRDS("hhv_france/clean_data/input_models/raw_data.RDS")
spatial_objects <-
  readRDS("hhv_france/clean_data/input_models/spatial_objects.RDS")$list_dpt


# Modifying label for sex --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  rename(sex = sexe) %>%
  mutate(sex = ifelse(sex == 'FEMININ', 'Female', 'Male'))




# Creating district variable ----------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  mutate(
    district = stringr::str_sub(dep_patient, start = 1, end = 2),
    district = ifelse(
      district == "97",
      stringr::str_sub(dep_patient, start = 1, end =
                         3),
      district
    ),
    city_sampling_site = dep_correspondant
  )

sort(as.numeric(unique(raw_data$district)))
district_excluded <- raw_data %>%
  filter(!district %in% spatial_objects$insee_dep)
unique(district_excluded$district)



# Renaming variable for IgG serology --------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  rename(
    results_hsv1 = hsvg_01,
    results_hsv2 = hsvg_04,
    results_cmv = cmvg_01,
    results_vzv = vzvg_01,
    results_ebv_eba = ebnae_01,
    results_ebv_vca = vcage_01
  )


# Label for test results-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  mutate_at(vars(starts_with("results_")),
            ~ case_when(. %in% c("p", "pf", "l", "n") ~ .,
                        TRUE ~ NA_character_))

# Removing patients with only NA results ----------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  filter(!(
    is.na(results_hsv1) &
      is.na(results_hsv2) &
      is.na(results_cmv) &
      is.na(results_vzv) &
      is.na(results_ebv_eba) &
      is.na(results_ebv_vca)
  ))



# Initial size of the dataset (for flowcharts)---------------------------------------------------------------------------------------------------------------------------------------------------------------------
fun_number_test <- function(df) {
  out <- df %>%
    select(
      one_of(
        "results_hsv1",
        "results_hsv2",
        "results_cmv",
        "results_vzv",
        "results_ebv_eba",
        "results_ebv_vca",
        "results_ebv"
      )
    ) %>%
    pivot_longer(
      cols = starts_with("results_"),
      names_to = "name",
      values_to = "values"
    ) %>%
    filter(!is.na(values)) %>%
    group_by(name) %>%
    count() %>%
    ungroup() %>%
    mutate(N = sum(n))
  
  print(out)
}

print("Initial dataset")
fun_number_test(raw_data)

# Filtering out individuals for which we don't have any information -------------------------------------------------------------------------------------------------------------------------------
file_filtered <- raw_data %>%
  filter(district %in% spatial_objects$insee_dep)

print("FILTERING OUT SAMPLING WITHOUT GEOGRAPHIC INFORMATION")
print("FILTERING OUT SAMPLING FROM PATIENTS LIVING ABROAD")
fun_number_test(raw_data %>%
                  filter(!district %in% spatial_objects$insee_dep))

table(file_filtered$district)



# Removing strange age ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("FILTERING OUT SAMPLING FROM PATIENT WITH NEGATIVE OR MISSING AGE")
fun_number_test(file_filtered %>%
                  mutate(age = as.numeric(age)) %>%
                  filter(age < 0 | is.na(age)))


file_filtered <- file_filtered %>%
  mutate(age = as.numeric(age)) %>%
  filter(age >= 0 & !is.na(age))

min(file_filtered$age)
max(file_filtered$age)

# Creating column year ----------------------------------------------------------------------------------------------------------------------------------------------------
file_filtered$date_sampling <-
  as.Date(file_filtered$date_prelevement, tryFormats = c("%d/%m/%Y"))

file_filtered %>% select(date_sampling, date_prelevement) %>%  view


min(file_filtered$date_sampling, na.rm = T)
max(file_filtered$date_sampling, na.rm = T)


# Filtering over study period ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("FILTERING OUT SAMPLING WITH MISSING DATE OF SAMPLING")
fun_number_test(file_filtered %>%
                  filter(is.na(date_sampling)))
print("FILTERING OUT SAMPLING OUTSIDE OF THE STUDY PERIOD")
fun_number_test(file_filtered %>%
                  filter(date_sampling < start_obs |
                           date_sampling > "2022-12-31"))

file_filtered <- file_filtered %>%
  filter(!is.na(date_sampling)) %>%
  filter(date_sampling >= start_obs & date_sampling <= "2022-12-31")

# Columns selection ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_filtered <- file_filtered %>%
  dplyr::select(-colonne1) %>%
  dplyr::select(
    num_patient,
    age,
    sex,
    district,
    city_sampling_site,
    date_prelevement,
    date_sampling,
    results_hsv1,
    results_hsv2,
    results_cmv,
    results_vzv,
    results_ebv_eba,
    results_ebv_vca
  )


# Adding time index -------------------------------------------------------------------------------------------------------------------------------------------------------------------

id_time <- seq(as.Date(start_obs),
               as.Date("2022-12-31"),
               by = "days") %>%
  tibble() %>%
  mutate(id_time = row_number())


colnames(id_time)[1] <- "date_sampling"

file_filtered <- file_filtered %>%
  left_join(id_time, by = "date_sampling")



# Adding district name ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects$geometry <- NULL

file_filtered <- file_filtered %>%
  left_join(spatial_objects %>%
              rename(district = insee_dep), by = "district")



# Individual numeric ID -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
id_patient_df <- file_filtered %>%
  distinct(num_patient) %>%
  mutate(id_patient = row_number())


file_filtered <- file_filtered %>%
  left_join(id_patient_df,
            by = "num_patient")

rm("id_patient_df")

file_filtered <- file_filtered %>%
  dplyr::select(-num_patient) %>%
  dplyr::relocate(id_patient, .before = age)



# Creating EBV outcome ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_filtered <- file_filtered %>%
  mutate(
    results_ebv = case_when(
      str_detect(results_ebv_eba,"p") ~ "p",
      str_detect(results_ebv_vca,"p") ~ "p",
      str_detect(results_ebv_eba,"n") & str_detect(results_ebv_vca,"n")~ "n",
      str_detect(results_ebv_eba,"n") &  str_detect(results_ebv_vca,"l")  ~ "n",
      str_detect(results_ebv_vca,"n") & str_detect(results_ebv_eba,"l")  ~ "n",
      str_detect(results_ebv_eba,"n") &  is.na(results_ebv_vca)  ~ "n",
      str_detect(results_ebv_vca,"n") & is.na(results_ebv_eba)  ~ "n",
      str_detect(results_ebv_vca,"l") & is.na(results_ebv_eba)  ~ "l",
      str_detect(results_ebv_eba,"l") & is.na(results_ebv_vca)  ~ "l",
      TRUE ~ NA_character_
    )
  )

#Checking it's fine
file_filtered %>%
  group_by(results_ebv_eba,
           results_ebv_vca,
           results_ebv) %>%
  count() 

#Yes it is


# Pivoting data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fun_number_test(file_filtered)


#Removng useless rows
file_filtered <- file_filtered %>%
  pivot_longer(cols = c(results_hsv1,
                        results_hsv2,
                        results_vzv,
                        results_cmv,
                        results_ebv),
               names_to = "name",
               values_to = "value") 


# Removing useless rows---------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_filtered <- file_filtered %>%
  filter(!is.na(value))


# Removing inconclusive test ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_filtered <- file_filtered %>% 
filter(value != "l") 
  
##Data when not keeping only 1 test/patient - used for sensitivity analysis
saveRDS(file_filtered %>%
          mutate(id_row_init = row_number()),
        "hhv_france/clean_data/input_models/data_clean_for_betabinomial.RDS")

# Keeping one test/ind-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Keeping only 1 test")
file_filtered %>%
  group_by(name, id_patient) %>%
  arrange(id_time) %>%
  filter(row_number() != 1) %>%
  ungroup() %>%
  group_by(name) %>%
  count() %>%
  ungroup() %>%
  filter(name != "results_ebv") %>%
  mutate(total = sum(n))



file_filtered <- file_filtered %>%
  group_by(name, id_patient) %>%
  arrange(id_time) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(id_row_init = row_number())



print("COMPARING INITIAL AND FINAL DATASET")
fun_number_test(raw_data)

file_filtered %>%
  group_by(name) %>%
  count()

file_filtered %>%
  group_by(name) %>%
  count() %>%
  ungroup() %>%
  mutate(id_row_init = row_number()) %>%
  mutate(total = sum(n)) %>%
  select(-id_row_init)

file_filtered %>% 
  filter(name == "results_ebv") %>% 
  filter(!is.na(results_ebv_vca)) %>% nrow
file_filtered %>% 
  filter(name == "results_ebv") %>% 
  filter(!is.na(results_ebv_eba)) %>% nrow

file_filtered %>% 
  filter(name != "results_ebv") %>% 
  dplyr::select(name,
                value) %>% 
  rbind(file_filtered %>% 
          filter(name == "results_ebv") %>% 
          filter(!is.na(results_ebv_vca)) %>% 
          dplyr::transmute(name = 'results_ebv_vca',
                           value = results_ebv_vca)) %>% 
  rbind(file_filtered %>% 
          filter(name == "results_ebv") %>% 
          filter(!is.na(results_ebv_eba)) %>% 
          dplyr::transmute(name = 'results_ebv_eba',
                           value = results_ebv_vca)) %>% 
  .$name %>% length()

# Saving ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(file_filtered,
        "hhv_france/clean_data/input_models/data_clean.RDS")
