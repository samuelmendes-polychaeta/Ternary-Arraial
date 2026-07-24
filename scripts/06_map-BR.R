library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(geobr)

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
