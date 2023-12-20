
library(tidyverse)
library(sf)

shape_river <- st_read("ecrins/ecrins_network.shp")
shape_basin <- st_read("Basin/Export_Output_11.shp")
shape_dams <- st_read("Barriers/Hampshire_Avon.shp")

source("functions_network_creation.R")