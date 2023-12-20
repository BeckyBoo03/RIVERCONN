# Case Study 2.2 - Bristol Avon Basin ####

# The first case study using RIVERCONN

# 1. Introduction ####

# 1.1 Packages needed ####

# Make sure those libraries are installed, updated, and loaded

# Install the latest version of RIVERCONN
# library(devtools)
# install_github("damianobaldan/riverconn")

library(tidyverse)
library(sf)
library(raster)
library(ggspatial)
library(viridis)
library(igraph)
library(riverconn)
library(elevatr)
library(gridExtra)
library(ggnetwork)
library(lwgeom)
library(gridExtra)
library(corrmorant)
library(RANN)
library(ggpubr)
library(cowplot)
library(extrafont)

# 1.3 Sources of the spatial data ####

# Load the data used for the analysis
shape_river<- st_read("Case Study 2.2/Bristol_Avon_Basin/River_Network/Export_Output_3.shp")
shape_basin <- st_read("Case Study 2.2/Bristol_Avon_Basin/Basin/Export_Output_2.shp")
shape_dams <- st_read("Case Study 2.2/Bristol_Avon_Basin/Barriers/Export_Output_4.shp")

# 1.4 Functions to support data processing ####

# This file includes the functions created for the RIVERCONN code
source("functions_network_creation.R")

# 2. Data Pre-processing ####

# 2.1 Shape files processing ####

# Set threshold for the length of river reach included in the analysis
threshold = 20

# Prune HydroRIVERS network
shape_river_small <- shape_river[as.numeric(shape_river$UPLAND_SKM) > threshold,]

# The dams shapefile can be converted to a shapefile and a unique id can be given 
# to each point
dams_to_points <- shape_dams %>%
  st_as_sf %>%
  st_centroid %>%
  mutate(id = 1:nrow(.))%>% 
  dplyr::select(id)%>% 
  st_transform(crs = "+proj=longlat +datum=WGS84")

# The shape files can be plot to check there are no errors in the geolocalization 
# of the data
ggplot() +
  coord_fixed() +
  theme_minimal() +
  ggspatial::layer_spatial(shape_basin, fill = NA, 
                           color = "gray25") +
  ggspatial::layer_spatial(shape_river_small, 
                           aes(color = log10(UPLAND_SKM)) )+
  ggspatial::layer_spatial(dams_to_points) +
  scale_color_viridis(direction = -1, 
                      name= "upstream area \n(log10[Km^2])",
                      limits = c(1, 4),
                      breaks = c(0, 0.5,1,1.5,2,2.5,3,3.5,4)) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.key.height = unit(2, "cm")) +
  ggspatial::annotation_scale(location = "bl", 
                              style = "ticks") +
  ggspatial::annotation_north_arrow(location = "br")+
  labs(caption = "Black dots are the position of the dams")+
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N"))

# 2.2 Confluences processing ####

# Simplify river shapefile
shape_river_simple <- shape_river_small %>%
  st_as_sf %>%
  st_union()

# Convert shapefile to point object
river_to_points <- shape_river_simple %>%
  st_as_sf %>%
  st_cast("POINT") %>%
  mutate(id = 1:nrow(.))

# Check points that overlap
joins_selection <- river_to_points %>%
  st_equals() %>%
  Filter(function(x){length(x) > 2}, .) %>%
  lapply(., FUN = min) %>%
  unlist() %>%
  unique()

# Filter original point shapefile to retain only confluences
river_joins <- river_to_points %>% 
  filter(id %in% joins_selection)

# Split polyline
shape_river_simplified <- lwgeom::st_split(shape_river_simple, river_joins) %>%
  st_collection_extract(.,"LINESTRING") %>%
  data.frame(id = 1:length(.), geometry = .) %>%
  st_as_sf() %>%
  mutate(length = st_length(.))

# Results can be checked with a plot.
ggplot() +
  coord_fixed() +
  ggspatial::layer_spatial(shape_river_simplified, 
                           aes(color = id))+
  scale_color_viridis(direction = -1, 
                      name= "Reach ID", 
                      limits = c(1, 69), 
                      breaks=c(1,10,20,30,40,50,60,69)) +
  ggspatial::layer_spatial(river_joins, shape = 1)+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(2, "cm"))+
  ggspatial::annotation_scale(location = "bl", 
                              style = "ticks") +
  ggspatial::annotation_north_arrow(location = "br")+
  labs(caption = "Hollow points are the position of the junctions")+
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N"))

# 2.3 Dams processing ####

# Snap dams
dams_snapped <- snap_to_river(dams_to_points,
                              shape_river_simple %>% st_sf(),
                              max_dist = 1000)

# Retain dams that were snapped
dams_snapped_reduced <-dams_snapped[st_contains(shape_river_simple %>% st_sf(), 
                                                dams_snapped,
                                                prepared = FALSE)[[1]],]

# Check if dams were snapped properly to the network (all the distances should 
# be zero)
st_distance(dams_snapped_reduced, 
            shape_river_simple) %>% sum

dams_snapped_reduced_joined <- dams_snapped_reduced %>%
  mutate(cluster =
           st_equals(dams_snapped_reduced, dams_snapped_reduced) %>%
           sapply(., FUN = function(x) paste0(x, collapse = "_"))) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(id_dam = as.character(1:nrow(.)), pass_u = 0.4, pass_d = 0.8) %>%
  as.data.frame %>%
  st_as_sf() %>%
  st_join(., shape_river_small, join = st_is_within_distance, dist = 10 ) %>% 
  group_by(id) %>%
  slice(which.max(UPLAND_SKM)) %>% 
  ungroup()

headwaters_checking <- headwaters_dam(dams_snapped_reduced_joined, 
                                      shape_river_simple)
head(headwaters_checking$flag_headwater)

# Check with plot the position of dams
my_plot <- ggplot() +
  coord_fixed() +
  ggspatial::layer_spatial(shape_river_simplified, color = "#6699FF")+
  ggspatial::layer_spatial(dams_snapped_reduced, shape = 16) +
  ggspatial::layer_spatial(shape_basin, fill = NA, color = "gray70") +
  theme_minimal() +
  ggspatial::annotation_north_arrow(location = "br") +
  ggspatial::annotation_scale(location = "bl", style = "ticks") +
  labs(caption = "Black Dots are the Position of Barriers") +
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N")) +
  theme(text = element_text(size = 20))

# Save the plot as a PNG image with a transparent background
ggsave("my_ggplot.png", plot = my_plot, width = 8, height = 6, bg = "transparent")


# 2.4 Shape file slicing ####

# Create junction point shapefile
network_links <- rbind(
  dams_snapped_reduced_joined %>% 
    mutate(type = "dam", id_barrier = id_dam) %>%
    dplyr::select(type, id_barrier, pass_u, pass_d),
  river_joins %>% mutate(type = "joint") %>%
    dplyr::select(type) %>%
    mutate(id_barrier = NA, pass_u = NA, pass_d = NA) %>%
    rename(geometry = x)) %>%
  mutate(id_links = 1:nrow(.))

# Split river network
river_net_simplified <- lwgeom::st_split(shape_river_simple, network_links) %>%
  st_collection_extract(.,"LINESTRING") %>%
  data.frame(NodeID = 1:length(.), geometry = .) %>%
  st_as_sf() %>%
  mutate(length = st_length(.)) %>%
  st_join(., shape_river_small, join = st_is_within_distance, dist = 0.01 ) %>% 
  group_by(NodeID) %>%
  slice(which.max(UPLAND_SKM)) %>% 
  ungroup()

# A plot can reveal if the process was successful.
ggplot() +
  coord_fixed() +
  ggspatial::layer_spatial(river_net_simplified, 
                           color = "gray70")+
  ggspatial::layer_spatial(network_links, 
                           aes(shape = type))+
  scale_shape_manual(name = "Splitting points", 
                     values=c("dam" =17,"joint" = 23))+
  theme_minimal() +
  theme(legend.position = "bottom")+
  ggspatial::annotation_north_arrow(location = "br")+
  ggspatial::annotation_scale(location = "bl", style = "ticks")+
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N"))


library(ggplot2)

# Your ggplot code
my_plot <- ggplot() +
  coord_fixed() +
  ggspatial::layer_spatial(river_net_simplified, 
                           color = "gray70")+
  ggspatial::layer_spatial(network_links, 
                           aes(shape = type))+
  scale_shape_manual(name = "Splitting points", 
                     values=c("dam" =17,"joint" = 23))+
  theme_minimal() +
  theme(legend.position = "bottom")+
  ggspatial::annotation_north_arrow(location = "br")+
  ggspatial::annotation_scale(location = "bl", style = "ticks")+
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N"))

# Save the plot as a PNG image with a transparent background
ggsave("my_ggplot.png", plot = my_plot, width = 8, height = 6, bg = "transparent")



confluences <- multiple_confluences(river_net_simplified)
head(confluences)

shp_check <- check_components(network_links, river_net_simplified)
head(shp_check)

# 2.5 Adding additional attributes to the river shapefile ####

# get DEM and transform to data frame with coordinates
elevation <- get_elev_raster(shape_basin, z = 8)
catchment_DEM <- raster::as.data.frame(elevation, xy = TRUE)

# Get coordinates of the river network segments
river_net_simplified_centroids <- river_net_simplified %>%
  st_as_sf() %>%
  st_centroid()

# Get coordinates of both elements for joining
Coord_Edges <- st_coordinates(river_net_simplified_centroids)
Coord_DEM <- catchment_DEM[,1:2]

# Matching each centroid with its closer altitude point to later obtain the 
# altitudes
matching_altitudes <- RANN::nn2(data=Coord_DEM, 
                                query = Coord_Edges, 
                                k=1, radius = 1)[[1]] 

# Get values and add to the river shapefile
catchment_DEM <- catchment_DEM[matching_altitudes,3]
river_net_simplified <- river_net_simplified %>% mutate(alt = catchment_DEM)

# To avoid issues in the network creation process, retain only the important 
# columns
river_net_simplified <- river_net_simplified %>% 
  dplyr::select(NodeID, length, alt, DIST_DN_KM, UPLAND_SKM)

# Plotting can show weird behaviours of the join operation.
my_plot <- ggplot() +    
  coord_fixed() +
  ggspatial::layer_spatial(river_net_simplified, 
                           aes(color = alt))+
  scale_color_viridis(name = "Elevation")+
  theme_minimal() +
  theme(legend.position = "bottom")+
  ggspatial::annotation_north_arrow(location = "br")+
  ggspatial::annotation_scale(location = "bl", style = "ticks")+
  scale_x_continuous(labels = function(x) paste0(x, "\u00B0 W")) +
  scale_y_continuous(labels = function(y) paste0(y, "\u00B0 N"))
# Save the plot as a PNG image with a transparent background
ggsave("elevation.png", plot = my_plot, width = 8, height = 6, bg = "transparent")


# 3. Creation of the igraph object ####

outlet <- river_net_simplified$NodeID[river_net_simplified$DIST_DN_KM == 0 ]

river_graph <- create_network(network_links, river_net_simplified, outlet)

# Check igraph object
river_graph

# check river_graph edges
igraph::get.edge.attribute(river_graph) %>% names

# check river_graph vertices
igraph::get.vertex.attribute(river_graph) %>% names

# 3.1 Editing the igraph object ####

# update length attribute
V(river_graph)$length <- V(river_graph)$length / 10000
hist(V(river_graph)$length)

# update length attribute
V(river_graph)$name <- as.character(V(river_graph)$name)

# Function for organism that prefers high elevation
suit_fun_high <- function(x){dnorm(x, mean = 1500, sd = 500)*410*sqrt(3*pi)}

# Function for organism that prefers low elevation
suit_fun_low <- function(x){exp(- 0.001*x)}

# Calculate HSI for the igraph nodes
V(river_graph)$HSI_low <- suit_fun_low(V(river_graph)$alt)
V(river_graph)$HSI_high <- suit_fun_high(V(river_graph)$alt)

# Calculate weighted usable length for igraph nodes
V(river_graph)$WUL_low <- V(river_graph)$HSI_low * V(river_graph)$length
V(river_graph)$WUL_high <- V(river_graph)$HSI_high * V(river_graph)$length

# Plot the two response functions
1:1:220 %>%
  data.frame("Elevation" = .,
             "Low" = suit_fun_low(.), 
             "High" = suit_fun_high(.)) %>%
  pivot_longer(c("Low", "High"), 
               names_to = "Type", values_to = "HSI") %>%
  ggplot() + 
  geom_line(aes(x = Elevation, y = HSI, color = Type))+
  theme_bw() + xlab("Elevation (m)") + ylab("Habitat Suitability Index (HSI)")

# 3.2 Plotting the igraph object ####

# The centroids of the reaches are extracted. The habitat suitabilities are 
# plotted.

# Extract reaches centroids
river_net_simplified_centroids <- river_net_simplified %>%
  st_as_sf() %>%
  st_centroid() 

# get the centroids coordinates
coordinates <- river_net_simplified_centroids %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,1], 
                lon = sf::st_coordinates(.)[,2]) %>%
  dplyr::select(lat, lon) %>%
  st_set_geometry( NULL)

# fortify the igraph object
gg0 <- ggnetwork(river_graph, 
                 layout = coordinates %>% as.matrix(), 
                 scale = FALSE)

grid.arrange(
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = HSI_high)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("High-altitude organism") +
    labs(caption = "Network directionality not shown"), 
  
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = WUL_high)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("High-altitude organism") +
    labs(caption = "Network directionality not shown"), 
  
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = HSI_low)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("Low-altitude organism") +
    labs(caption = "Network directionality not shown"),
  
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = WUL_low)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("Low-altitude organism") +
    labs(caption = "Network directionality not shown"),
  ncol=2, nrow=2)


my_plot <- ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_fixed() +
  geom_nodes(aes(color = HSI_high)) +
  geom_edges(alpha = 0.5) +
  scale_color_viridis()+
  theme_blank()+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.key = element_blank() #transparent legend key (background)
  )+
  labs(caption = "Network directionality not shown")
# Save the plot as a PNG image with a transparent background
ggsave("HSI.png", plot = my_plot, width = 8, height = 6, bg = "transparent")

my_plot <- ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_fixed() +
  geom_nodes(aes(color = HSI_low)) +
  geom_edges(alpha = 0.5) +
  scale_color_viridis()+
  theme_blank()+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.key = element_blank() #transparent legend key (background)
  )+
  labs(caption = "Network directionality not shown")
# Save the plot as a PNG image with a transparent background
ggsave("HSI.png", plot = my_plot, width = 8, height = 6, bg = "transparent")


# 4. Indices Calculations ####

# 4.1 Reach-scale Indices ####

# Initialize list where to store all the index calculation outputs
index <- list()
lab_index <- list()
letter_index <- list()

# 1: Symmetric Dendritic Connectivity Index (no biotic effects)
lab_index[[1]] <- "DCI with symmetric passabilities"
letter_index[[1]] <- "1"
index[[1]] <- index_calculation(graph = river_graph,                             # Input 'igraph' object
                               weight = "length",                                # graph vertex attribute to be used as w1 in equation 2.1 and 2.2
                               B_ij_flag = FALSE,                                # Defines if the contribution of Bij to CCI or RCI is calculated or not (eq. (2.3))
                               index_type = "reach",                             # Defines if CCI or RCI is calculated
                               index_mode = "from")                              # Defines if inbound or outbound links are used in RCI

# 2: Asymmetric Dendritic Connectivity Index (no biotic effects)
lab_index[[2]] <- "DCI with asymmetric passabilities"
letter_index[[2]] <- "2"
index[[2]] <- index_calculation(graph = river_graph,                             # Input 'igraph' object
                               weight = "length",                                # graph vertex attribute to be used as w1 in equation 2.1 and 2.2
                               B_ij_flag = FALSE,                                # Defines if the contribution of Bij to CCI or RCI is calculated or not (eq. (2.3))
                               dir_fragmentation_type = "asymmetric",            # Defines how equivalent barriers passability is calculated (eq (2.4) or 2.5)
                               index_type = "reach",                             # Defines if CCI or RCI is calculated
                               index_mode = "from")                              # Defines if inbound or outbound links are used in RCI

# Before calculating IIC a binary passability has to be defined
E(river_graph)$pass_u_bin <- ifelse(is.na(E(river_graph)$pass_u), NA, 0.4)       # I have changed this from the tutorial as the barriers can be passed in the case study
E(river_graph)$pass_d_bin <- ifelse(is.na(E(river_graph)$pass_d), NA, 0.8)       # I have changed this from the tutorial as the barriers can be passed in the case study

# 3: Symmetric Integral Index of Connectivity
lab_index[[3]] <- "IIC"
letter_index[[3]] <- "3"
index[[3]] <- index_calculation(graph = river_graph,                             # Input ‘igraph’ object
                               weight = "length",                                # Graph vertex attribute to be used as wi in eqs. (2.1) and (2.2)
                               pass_u = "pass_u_bin",                            # Defines upstream-to-downstream barriers passabilities 
                               pass_d = "pass_d_bin",                            # Defines downstream-to-upstream barriers passabilities 
                               param = 0.6,                                      # Defines the dispersal parameter PD in eq. (2.8) or 2.9 
                               index_type = "reach",                             # Defines if CCI or RCI is calculated
                               index_mode = "from")                              # Defines if inbound or outbound links are used in RCI

## Defining uniform weights for IIC
V(river_graph)$unif_w <- 1

# 4: Integral Index of Connectivity with uniform weights - only dams contribute to index value
lab_index[[4]] <- "IIC with uniform weights"
letter_index[[4]] <- "4"
index[[4]] <- index_calculation(graph = river_graph,                             # Input ‘igraph’ object
                               weight = "unif_w",                                # Only the presence of barriers contribute to river fragmentation
                               dir_fragmentation_type = "asymmetric",            # Defines how equivalent barriers passability is calculated (eq (2.4) or 2.5)
                               pass_u = "pass_u_bin",                            # Defines upstream-to-downstream barriers passabilities 
                               pass_d = "pass_d_bin",                            # Defines downstream-to-upstream barriers passabilities 
                               param = 0.6,                                      # Defines the dispersal parameter PD in eq. (2.8) or 2.9 
                               index_type = "reach",                             # Defines if CCI or RCI is calculated
                               index_mode = "from")                              # Defines if inbound or outbound links are used in RCI
 
# 5: Population Connectivity Index for lowland fish
lab_index[[5]] <- "PCI with symmetric dispersal"
letter_index[[5]] <- "5"
index[[5]] <- index_calculation(graph = river_graph,                             # Input ‘igraph’ object
                                weight = "WUL_low",                              # Graph vertex attribute to be used as wi in eqs. (2.1) and (2.2)
                                param = 0.6,                                     # Defines the dispersal parameter PD in eq. (2.8) or 2.9 
                                index_type = "reach",                            # Defines if CCI or RCI is calculated
                                index_mode = "from")                             # Defines if inbound or outbound links are used in RCI

# 6: Population Connectivity Index for lowland passive drifter
lab_index[[6]] <- "PCI with assymetric dispersal"
letter_index[[6]] <- "6"
index[[6]] <- index_calculation(graph = river_graph,                             # Input ‘igraph’ object
                                weight = "WUL_low",                              # Graph vertex attribute to be used as wi in eqs. (2.1) and (2.2)
                                dir_distance_type  = "asymmetric",               # Defines if directionality is accounted for in Bij calculation (eqs. (2.8) - (2.11))
                                param_u = 0.4,                                   # Defines the upstream dispersal parameter PDu in eq. (2.10) or 2.11 
                                param_d = 0.6,                                   # Defines the downstream dispersal parameter PDd
                                index_type = "reach",                            # Defines if CCI or RCI is calculated
                                index_mode = "from")                             # Defines if inbound or outbound links are used in RCI

# A quick check to the variation range
data.frame("Baseline CCI" = do.call(rbind, lab_index),
           "min" = t(sapply(index, sapply, min))[,4], # baseline in is column 4
           "mean" = t(sapply(index, sapply, mean))[,4],
           "max" = t(sapply(index, sapply, max))[,4] )

# Resulting maps can be plotted to visually compare the results. 

# Initialize empty list
plot_list <- list()

# iterate through list length
for (i in 1:length(index)) {
  
  # Join d_i information with dams shapefile
  river_plot <- river_net_simplified %>%
    mutate(name = as.character(NodeID)) %>%
    left_join(index[[i]])
  
  # plot
  plot_iter <- ggplot() +
    coord_fixed()+
    ggspatial::layer_spatial(river_plot, aes(color = index))+ 
    scale_color_viridis(name = paste0("Index Value"),
                        direction = -1, 
                        breaks = seq(min(index[[i]]$index), 
                                     max(index[[i]]$index), 
                                     length.out=5))+
    theme_void()+
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.key.height = unit(0.8, "cm"),
          legend.key.size = unit(0.3, "cm"),
          plot.margin = margin(5, 5, 5, 5),
          text = element_text(size = 13)) + 
    ggtitle(paste0(letter_index[[i]], ") ", lab_index[[i]])) 
  plot_list[[i]] <- plot_iter}
plot_final <- ggpubr::ggarrange(plotlist = plot_list, nrow = 3, ncol = 2, legend = "right")
plot_final
ggsave("plot_final.png", plot = plot_final, width = 10, height = 10, bg = "transparent")

# A correlation plot helps identifying the correlated indices
d_index_chart <- data.frame("d_i_1" = index[[1]]$index, 
                            "d_i_2" = index[[2]]$index,
                            "d_i_3" = index[[3]]$index,
                            "d_i_4" = index[[4]]$index,
                            "d_i_5" = index[[5]]$index,
                            "d_i_6" = index[[6]]$index)

colnames(d_index_chart) <- c(letter_index[[1]], 
                             letter_index[[2]], 
                             letter_index[[3]], 
                             letter_index[[4]], 
                             letter_index[[5]], 
                             letter_index[[6]])

# Spearman correlations between rankings calculated with different priorities
cor(d_index_chart, method = "spearman")

corrmorant::corrmorant(d_index_chart, 
                       corr_method = "spearman", 
                       style = "binned")+
  theme(legend.position = "bottom") +
  scale_color_viridis(name = "Correlation", 
                      option = "E", 
                      direction = -1, 
                      limits = c(-1, 1)) +
  theme(axis.text = element_text(size=5),
        axis.text.x = element_text(angle = 45) )


# 4.2 Barriers Prioritization ####

# Dams metadata is to be defined
barriers_metadata <- data.frame("id_barrier" =  E(river_graph)$id_barrier
                                [!is.na(E(river_graph)$id_barrier)],
                                "pass_u_updated" = 1,
                                "pass_d_updated" = 1)
head(barriers_metadata)

# Next, the dams prioritization can be executed. 
# Initialize list where to store all the index calculation outputs
d_index <- list()
lab_d_index <- list()
letter_d_index <- list()

# 1: Symmetric Dendritic Connectivity Index (no biotic effects)
lab_d_index[[1]] <- "DCI with symmetric passabilities"
letter_d_index[[1]] <- "1"
d_index[[1]] <- d_index_calculation(graph = river_graph, 
                                    barriers_metadata = barriers_metadata,
                                    B_ij_flag = FALSE, 
                                    parallel = FALSE)

# 2: Asymmetric Dendritic Connectivity Index (no biotic effects)
lab_d_index[[2]] <- "DCI with asymmetric passabilities"
letter_d_index[[2]] <- "2"
d_index[[2]] <- d_index_calculation(graph = river_graph, 
                                    barriers_metadata = barriers_metadata, 
                                    B_ij_flag = FALSE, 
                                    parallel = FALSE, 
                                    dir_fragmentation_type = "asymmetric")

## Before calculating IIC a binary passability has to be defined
E(river_graph)$pass_u_bin <- ifelse(is.na(E(river_graph)$pass_u), NA, 0.4)
E(river_graph)$pass_d_bin <- ifelse(is.na(E(river_graph)$pass_d), NA, 0.8)

# 3: Symmetric Integral Index of Connectivity
lab_d_index[[3]] <- "IIC"
letter_d_index[[3]] <- "3"
d_index[[3]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low",
                                    pass_u = "pass_u_bin",
                                    pass_d = "pass_d_bin",
                                    param = 0.6,
                                    parallel = FALSE)

## Defining uniform weights for IIC
V(river_graph)$unif_w <- 1

# 4: Integral Index of Connectivity with uniform weights
lab_d_index[[4]] <- "IIC with uniform weights"
letter_d_index[[4]] <- "4"
d_index[[4]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "unif_w",
                                    dir_fragmentation_type = "asymmetric",
                                    pass_u = "pass_u_bin",
                                    pass_d = "pass_d_bin",
                                    param = 0.6,
                                    parallel = FALSE)

# 5: Population Connectivity Index for lowland fish
lab_d_index[[5]] <- "PCI with symmetric dispersal"
letter_d_index[[5]] <- "5"
d_index[[5]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low", 
                                    param = 0.6,
                                    parallel = FALSE) 

# 6: Population Connectivity Index for lowland passive drifter
## note that units of param_d must be consistent with length field (10s of km)
lab_d_index[[6]] <- "PCI with asymmetric dispersal"
letter_d_index[[6]] <- "6"
d_index[[6]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low", 
                                    dir_distance_type  = "asymmetric", 
                                    param_u = 0.4,
                                    param_d = 0.6,
                                    parallel = FALSE)

# A quick check to the variation range of CCI for the baseline.
data.frame("Baseline CCI" = do.call(rbind, lab_d_index),
           "min" = t(sapply(d_index, sapply, min))[,4],
           "mean" = t(sapply(d_index, sapply, mean))[,4],
           "max" = t(sapply(d_index, sapply, max))[,4] )

# Initialize empty list
plot_list <- list()

# iterate through list length
for (i in 1:length(d_index)) {
  # Join d_i information with dams shapefile
  network_links_plot <- network_links %>%
    filter(type == "dam") %>% 
    left_join(d_index[[i]]) %>%
    mutate(d_rank = rank(desc(d_index)))
  # plot
  plot_iter <- ggplot() +
    coord_fixed() +
    ggspatial::layer_spatial(river_net_simplified, color = "gray70")+
    ggspatial::layer_spatial(network_links_plot,aes(color = d_rank, size = 1/d_rank), alpha = 0.8)+
    scale_color_viridis(name = "Ranking",direction = -1, limit= c(1, 154), breaks = c(1,50,100,154))+
    theme_void()+
    guides(size = FALSE) +
    theme(legend.direction = "vertical",
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(3, "cm"),
          text = element_text(size = 13))+
    ggtitle(paste0(letter_d_index[[i]], ") ", lab_d_index[[i]]))
  plot_list[[i]] <- plot_iter
}
plot_final <- ggpubr::ggarrange( plotlist = plot_list,common.legend = TRUE, legend = "right")
plot_final
ggsave("plot_final.png", plot = plot_final, width = 13, height = 7, bg = "transparent")


# A correlation plot helps identifying the correlated indices
d_index_chart <- data.frame("d_i_1" = d_index[[1]]$d_index, 
                            "d_i_2" = d_index[[2]]$d_index,
                            "d_i_3" = d_index[[3]]$d_index,
                            "d_i_4" = d_index[[4]]$d_index,
                            "d_i_5" = d_index[[5]]$d_index,
                            "d_i_6" = d_index[[6]]$d_index)

colnames(d_index_chart) <- c(letter_d_index[[1]], 
                             letter_d_index[[2]], 
                             letter_d_index[[3]], 
                             letter_d_index[[4]], 
                             letter_d_index[[5]], 
                             letter_d_index[[6]])

# Spearman correlations between rankings calculated with different priorities
cor(d_index_chart, method = "spearman")

corrmorant::corrmorant(d_index_chart, 
                       corr_method = "spearman", 
                       style = "binned")+
  theme(legend.position = "bottom") +
  scale_color_viridis(name = "Correlation", 
                      option = "E", 
                      direction = -1, 
                      limits = c(-1, 1)) +
  theme(axis.text = element_text(size=5),
        axis.text.x = element_text(angle = 45))