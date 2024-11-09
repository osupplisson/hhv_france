#Cleaning R objects ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(xlsx)
library(janitor)


# Importing data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_1 <- read_excel("hhv_france/raw_data/extraction_2/GLPI_30497_onyme1.xlsx",
                                            col_types = rep("text", 41)) %>%
  clean_names()

file_2 <- read_excel("hhv_france/raw_data/extraction_2/GLPI_30497_onyme2.xlsx",
                                            col_types = rep("text", 41)) %>%
  clean_names()

file_3 <- read_excel("hhv_france/raw_data/extraction_3/glpi_32764 V3_onyme.xlsx",
                     col_types = rep("text", 41)) %>%
  clean_names() %>% 
  rename(date_prelevement = `x12`)

file <- rbind(file_1,
              file_2,
              file_3)

table(format(as.Date(file_1$date_prelevement, tryFormats = c("%d/%m/%Y")), "%Y"))
table(format(as.Date(file_2$date_prelevement, tryFormats = c("%d/%m/%Y")), "%Y"))
table(format(as.Date(file_3$date_prelevement, tryFormats = c("%d/%m/%Y")), "%Y"))

# save --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(file, "hhv_france/clean_data/input_models/raw_data.RDS")
