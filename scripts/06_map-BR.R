library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(geobr)
library(here)
library(openxlsx)

#coords
coords<-read.xlsx(here("tables", "data.xlsx"), sheet = "coords") 

# Países da América do Sul
america_sul <- ne_countries(
  continent = "South America",
  returnclass = "sf"
)

# Brasil (apenas um polígono, sem estados)
brasil <- america_sul[america_sul$admin == "Brazil", ]

# Estado do Rio de Janeiro
rj <- read_state(code_state = "RJ", year = 2020)

ggplot() +
  
  # Países vizinhos
  geom_sf(data = america_sul,
          fill = "#C2AD95",
          color = "black",
          linewidth = 0.3) +
  
  # Brasil
  geom_sf(data = brasil,
          fill = "#C2AD95",
          color = "black",
          linewidth = 0.3) +
  
  # Contorno do Rio de Janeiro
  geom_sf(data = rj,
          fill = NA,
          color = "red",
          linewidth = 0.3) +
  
  # Barra de escala
  annotation_scale(
    location = "br",
    width_hint = 0.22
  ) +
  
  # Norte
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  
  coord_sf(
    xlim = c(-85, -30),
    ylim = c(-36, 13),
    expand = FALSE
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#9ED2DE"),
    panel.grid.major = element_line(color = "#E8EDED", linewidth = 0.3),
    panel.grid.minor = element_line(color = "black", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# 1. Load your coordinate data
coords <- read.xlsx(here("tables", "data.xlsx"), sheet = "coords")

# 2. Define the geographic bounding box (converted S and W to negative numbers)
lat_min <- -23.03
lat_max <- -22.94
lon_min <- -42.04
lon_max <- -41.96

# 3. Load high-detail coastline/land boundary
land <- ne_countries(scale = 10, returnclass = "sf")

# 4. Create the Map using ggplot2
ggplot() +
  # Land Polygon
  geom_sf(data = land, fill = "antiquewhite", color = "gray50", size = 0.4) +
  
  # Sampling Points
  geom_point(
    data = coords,
    aes(x = long, y = lat),
    color = "black", fill = "red", size = 3.5, shape = 21, stroke = 1
  ) +
  
  # Set Map Window
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  
  # Scale Bar & North Arrow
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    style = north_arrow_minimal()
  ) +
  
  # Plot Aesthetics
  labs(
    title = "Cabo Frio RESEX Sampling Points",
    x = "Longitude (°W)",
    y = "Latitude (°S)"
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "aliceblue"), # Ocean color
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

library(ggplot2)
library(sf)
library(geobr)
library(ggspatial)
library(openxlsx)
library(here)

# 1. Load spatial land data
rj <- read_state(code_state = "RJ", year = 2020)
brasil <- read_country(year = 2020)

# 2. Load your coordinates from Excel
# Assumes columns in Excel are named 'longitude' and 'latitude' (or 'lon'/'lat')
coords_df <- read.xlsx(here("tables", "data.xlsx"), sheet = "coords")

# Convert coordinate data frame to an sf object (WGS84 projection: EPSG 4326)
coords_sf <- st_as_sf(
  coords_df, 
  coords = c("long", "lat"), 
  crs = 4326
)

# 3. Build the plot
ggplot() +
  
  # Brazil landmass (covers the RJ region and surrounding land)
  geom_sf(
    data = brasil,
    fill = "#C2AD95",
    color = "black",
    linewidth = 0.3
  ) +
  
  # State boundary outline for Rio de Janeiro
  geom_sf(
    data = rj,
    fill = NA,
    color = NA,
    linewidth = 0.4
  ) +
  
  # Plot sampling points from coords dataset
  geom_sf(
    data = coords_sf,
    color = 'black',
    fill = "white",
    shape = 21,
    size = 5,
    stroke = 0.5
  ) +
  
  # Scale bar (bottom-right)
  annotation_scale(
    location = "br",
    width_hint = 0.25,
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm")
  ) +
  
  # North Arrow (top-right)
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.0, "cm"),
    width = unit(1.0, "cm"),
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm")
  ) +
  
  # Crop precisely to the requested Cabo Frio / RESEX bounding box
  coord_sf(
    xlim = c(-42.04, -41.96),
    ylim = c(-23.03, -22.94),
    expand = FALSE
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#9ED2DE"), # Ocean blue color
    panel.grid.major = element_line(color = "#E8EDED", linewidth = 0.3),
    panel.grid.minor = element_line(color = "black", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
