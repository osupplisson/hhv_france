#This code defines the population strata used during the post-stratifications step
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hhv_france/codes/00_packages.R")
library("readxl")


# Importing spatial data --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects <- readRDS("hhv_france/clean_data/input_models/spatial_objects.RDS")$list_dpt


# Function for importing pop ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
importing_population <- function(year_choice) {
  df <- read_excel("global_raw_data/data_departements/estim-pop-dep-sexe-aq-1975-2023.xls",
    sheet = year_choice, skip = 4
  ) %>%
    rename_all(~ str_replace(., "\\.\\.\\.", "\\-")) %>%
    as_tibble()

  # Changing names for 2 first cols
  colnames(df)[1] <- "insee_dep"
  colnames(df)[2] <- "nom"
  # Dropping na
  df <- df %>% filter(!is.na(nom))

  # Merging Corsica
  df <- df %>%
    pivot_longer(
      cols = -c(
        insee_dep,
        nom
      ),
      names_to = "name",
      values_to = "population"
    ) %>%
    mutate(
      insee_dep = ifelse(insee_dep %in% c("2A", "2B"),
        "20",
        insee_dep
      ),
      nom = ifelse(insee_dep == "20",
        "Corse",
        nom
      )
    ) %>%
    summarise(population = sum(population), .by = c(insee_dep, nom, name))

  # Renaming and creating sex cols
  df <- df %>%
    tidyr::separate_wider_delim(name,
      delim = "-",
      names = c("name", "order")
    ) %>%
    mutate(order = as.numeric(order)) %>%
    filter(name != "Total") %>%
    group_by(
      insee_dep,
      nom,
      name
    ) %>%
    arrange(order) %>%
    mutate(order_bis = row_number()) %>%
    mutate(
      sex = case_when(
        order_bis == 1 ~ "All",
        order_bis == 2 ~ "Male",
        order_bis == 3 ~ "Female"
      )
    )%>%
    ungroup() %>% 
    mutate(year = year_choice)

  return(df)
}




# Applying function -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
year_set <- as.character(seq(2018,2022,1))
for(year_choice in year_set){
  for.loop <-importing_population(year_choice) 
  if(year_choice == year_set[1]){
   df_year <- for.loop 
  }else{
    df_year <- rbind(for.loop, 
          df_year)
  }
}


# Age ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_year <- df_year %>% 
  dplyr::select(-order_bis) %>% 
  mutate(age_class = str_replace(name, " à ", "\\-"),
         age_class = str_remove(age_class, "ans"),
         age_class = str_replace(age_class, "  et plus" , "\\-\\+")) %>%
  tidyr::separate_wider_delim(age_class,
                              delim = "-",
                              names = c("lb", "ub")
  ) %>% 
  mutate(lb = as.numeric(lb),
         ub = as.numeric(ub)) 


df_year <- df_year %>% 
  dplyr::select(
    insee_dep,
    nom,
    sex,
    year,
    lb,
    ub,
    population
  )

df_year <- df_year %>% 
  mutate(
    age_class = case_when(
      lb == 0 ~ paste0(lb, ",", ub),
      !lb %in% c(0, 95) ~ paste0(lb - 1, ",", ub),
      lb == 95 ~ "95+"
    ),
    age_class = case_when(
      lb == 0 ~ paste0("[", age_class, "]"),
      !lb %in% c(0, 95) ~ paste0("(", age_class, "]"),
      lb == 95 ~ age_class
    )
  ) 

# Sex ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_year <- df_year %>% 
  filter(sex != "All")%>%
  mutate(sex_numeric = ifelse(sex == "Male", 1, 2)) %>% 
  rename(N = population)


# Adding long/lat ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects$geometry <- NULL

df_year <-df_year %>% 
  left_join(spatial_objects %>% 
              dplyr::select(insee_dep,
                            centroid_long,
                            centroid_lat,
                            district_numeric),
            by = c("insee_dep")) %>% 
  rename(id_district = insee_dep)



# Adding id vector --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
id_year <- tibble(year = unique(df_year$year)) %>% 
  arrange(year) %>% 
  mutate(year_numeric = row_number(),
         year_factor = factor(year_numeric, ordered = T)) 

id_age <- tibble(lb = unique(df_year$lb)) %>% 
  arrange(lb) %>% 
  mutate(age_numeric = row_number(),
         age_factor = factor(age_numeric, ordered = T)) 

df_strata <- df_year %>% 
  left_join(
    id_year, 
    by = "year"
  ) %>% 
  left_join(
    id_age,
    by = "lb"
  )


# Order -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_strata <- df_strata  %>% 
  mutate(sex_factor = factor(sex, levels = c("Male",
                                             "Female"))) %>% 
  dplyr::select(
    nom,
    centroid_long,
    centroid_lat,
    id_district,
    district_numeric,
    sex,
    sex_numeric,
    sex_factor,
    year,
    year_numeric,
    year_factor,
    lb,
    ub,
    age_class,
    age_numeric,
    age_factor,
    N
  ) 

df_strata %>% 
  group_by(sex,
           year) %>% 
  summarise(N = sum(N))

df_strata %>% 
  group_by(year) %>% 
  summarise(N = sum(N))

# Save --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(df_strata,
        "hhv_france/clean_data/input_models/strata.RDS")
