#This code produces the adjacency matrices further used for the BYM2 structured prior

# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(spdep)
library(janitor)
library(tidyverse)
library(tidylog)
library(dbscan)
# Importing shapefile from France -----------------------------------------------------------------------------------------------------------------------------------------------------------------
carto_dpt <-
  read_sf(
    "global_raw_data/data_for_cartography/1_DONNEES_LIVRAISON_2022-08-30/ADE_3-1_SHP_WGS84G_FRA/DEPARTEMENT.shp"
  ) %>%
  clean_names()



# Extracting Corsica ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
corsica <- carto_dpt %>%
  filter(insee_dep %in% c("2A", "2B"))

#Union of Corsican distritcs
union_corsica <- st_union(corsica$geometry[1],
                          corsica$geometry[2])

new_corsica <- corsica %>%
  filter(row_number() == 1) %>%
  mutate(
    id = 'XXX',
    nom_m = "CORSE",
    nom = 'Corse',
    insee_dep = "20",
    insee_reg = "94",
    geometry = union_corsica
  )


carto_dpt <- carto_dpt %>%
  rbind(new_corsica) %>%
  filter(!insee_dep %in% c("2A", "2B"))


rm("corsica")
rm("new_corsica")



# Continuous spatial basis ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extracting centroids
carto_dpt_centroids <- carto_dpt %>%
  mutate(centroids = st_centroid(st_geometry(.)))

# Geom to coordinate
carto_dpt_centroids <- carto_dpt_centroids %>%
  mutate(centroids.var = st_coordinates(centroids))



# Longitude basis ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
longitudes <- carto_dpt_centroids$centroids.var[, 1]
# Latitude basis ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
latitudes <- carto_dpt_centroids$centroids.var[, 2]

rm("carto_dpt_centroids")

# Binding -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_basis_dpt <- cbind(carto_dpt,
                           tibble("centroid_long" = longitudes),
                           tibble("centroid_lat" = latitudes))


# Neighboorhood matrix for CAR models---------------------------------------------------------------------------------------------------------------
# https://r-spatial.github.io/spdep/articles/nb.html
# Creating queen contiguity
queen_contiguity_nbobject <- poly2nb(as(carto_dpt, "Spatial"),
                                     queen = TRUE)

names(queen_contiguity_nbobject) <- carto_dpt$nom

# SPDED to INLA graph
spdep::nb2INLA(
  "hhv_france/clean_data/input_models/inla_graph_districts_queen.graph",
  queen_contiguity_nbobject
)


# Creating distance-based
long_lat_matrix_init <- spatial_basis_dpt %>%
  dplyr::select(nom,
                centroid_long,
                centroid_lat)

long_lat_matrix_init$geometry <- NULL
long_lat_matrix <- long_lat_matrix_init  %>%
  dplyr::select(centroid_long,
                centroid_lat) %>%
  as.matrix()


overseas <- long_lat_matrix_init %>%
  mutate(id = row_number()) %>%
  filter(nom %in% c("Guadeloupe",
                    "Martinique",
                    "Mayotte",
                    "La Réunion",
                    "Guyane"))

removing_overseas <- function(list,
                              overseas_neighbors = T) {
  list <- list %>%
    imap(function(.x, .y) {
      if (!.y %in% overseas$id) {
        .x[!.x %in% c(overseas$id)]
      } else if (.y %in% overseas$id) {
        .x[.x %in% c(overseas$id)]
      }
    })
  
  if (overseas_neighbors == T) {
    list <- list %>%
      imap(function(.x, .y) {
        if (.y %in% overseas$id) {
          c(overseas$id)[-which(c(overseas$id)==.y)]
        }
        else{
          .x
        }
      })
  }
  class(list) <- c("nb")
  list
}


# Delaunay triangulation neighbours and SOI neighbours are symmetric by design
delaunay_nbobject <- tri2nb(long_lat_matrix)

delaunay_nbobject <- delaunay_nbobject %>%
  removing_overseas()

spdep::nb2INLA(
  "hhv_france/clean_data/input_models/inla_graph_districts_delaunay.graph",
  delaunay_nbobject
)

soi_nbobject <-  graph2nb(soi.graph(tri2nb(long_lat_matrix), long_lat_matrix)) %>%
  removing_overseas()

spdep::nb2INLA(
  "hhv_france/clean_data/input_models/inla_graph_districts_soi.graph",
  soi_nbobject
)


#Gabriel
gab_nbobject <-
  graph2nb(gabrielneigh(long_lat_matrix), sym = TRUE) %>%
  removing_overseas()
spdep::nb2INLA(
  "hhv_france/clean_data/input_models/inla_graph_districts_gabriel.graph",
  gab_nbobject
)

#Relative graph
relative_nbobject <-
  graph2nb(relativeneigh(long_lat_matrix), sym = TRUE) %>%
  removing_overseas()
spdep::nb2INLA(
  "hhv_france/clean_data/input_models/inla_graph_districts_relative.graph",
  relative_nbobject
)


# Neighbors ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nb_list <- tibble(knn = seq(1, 40)) %>%
  split(as.factor(.$knn)) %>%
  imap(~knn2nb(knearneigh(long_lat_matrix, k = .x), row.names = spatial_basis_dpt$id)) %>%
  imap(~make.sym.nb(.x)) %>%
  imap(~removing_overseas(.x,overseas_neighbors = F))

names(nb_list) <- paste0("nb", names(nb_list))

nb_list %>%
  imap( ~ spdep::nb2INLA(
    paste0(
      "hhv_france/clean_data/input_models/inla_graph_districts_",
      .y,
      ".graph"
    ),
    .x
  ))


# ID numeric district -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_basis_dpt <- spatial_basis_dpt %>%
  mutate(district_numeric = row_number())



# Adding region names -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
region_name <-
  read_sf(
    "global_raw_data/data_for_cartography/1_DONNEES_LIVRAISON_2022-08-30/ADE_3-1_SHP_WGS84G_FRA/REGION.shp"
  ) %>%
  clean_names() %>%
  dplyr::select(nom_region = nom,
                insee_reg = insee_reg)
region_name$geometry <- NULL

spatial_basis_dpt <- spatial_basis_dpt %>%
  left_join(region_name, by = "insee_reg")



# Saving -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects <- list(
  list_dpt = spatial_basis_dpt,
  nb = list(
    'queen' = queen_contiguity_nbobject,
    'delaunay' = delaunay_nbobject,
    'soi' = soi_nbobject,
    'gabriel' = gab_nbobject,
    'relative' = relative_nbobject,
    'nb' = nb_list
  )
)

saveRDS(spatial_objects,
        "hhv_france/clean_data/input_models/spatial_objects.RDS")
