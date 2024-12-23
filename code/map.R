# Draw a map of the study area

library(raster)
library(ggspatial)
library(grid)
library(gridExtra)
library(magick)

source("code/model.R")

# main map
lims <- st_bbox(roi)

p_map <- 
  ggplot() +
  geom_sf(data = roads,
          colour = "darkgrey") +
  geom_sf(data = reservoir,
          fill = "skyblue") +
  geom_sf(data = roi,
          fill = "springgreen3") +
  geom_sf(data = coord,
          pch = 15,
          size = 0.1,
          colour = "black") +
  coord_sf(xlim = lims[c(1, 3)],
           ylim = lims[c(2, 4)],
           expand = FALSE) +
  annotation_scale(
    width_hint = 0.3,
    location = "br",
    bar_cols = c("grey60", "white"),
    style = "ticks"
  ) +
  annotation_north_arrow(
    location = "br", which_north = "true",
    pad_y = unit(0.4, "in")
  ) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "snow2"),
        panel.grid = element_blank())

# inset
sg <- 
  st_read("data/GIS/NationalMapPolygon.kml") %>% 
  # filter(FOLDERPATH == "Layers/Coastal_Outlines")
  filter(Name == "SINGAPORE - MAIN ISLAND")

study_box <- st_make_grid(roi, n = 1)

p_inset <- 
  ggplot() +
  geom_sf(data = sg, 
          colour = "black",
          fill = "snow2") +
  geom_sf(data = study_box,
          linejoin = "mitre", lineend = "butt",
          color = "red",
          fill = NA,
          lwd = 1) +
  coord_sf() +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_rect(fill = "white"))

# combine
pdf("fig/map.pdf", width = 6, height = 5)

grid.newpage()
mainmap <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) 
insetmap <- viewport(width = 0.3, height = 0.3, x = 0.3, y = 0.85)
print(p_map, vp = mainmap) 
print(p_inset, vp = insetmap)

dev.off()


# map plus sampling design
img1 <- image_read_pdf("fig/map.pdf")
# img2 <- image_read("fig/adaptive_cluster.png")
img2 <- image_read_pdf("fig/grid.pdf")

# Add labels to the images
img1_labeled <- 
  img1 %>% 
  image_resize("x800") %>% 
  image_annotate("A", size = 12, gravity = "northwest") 
img2_labeled <- 
  img2 %>% 
  image_resize("x400") %>% 
  image_border(color = "white", geometry = "0x50") %>%
  image_annotate("B", size = 12, gravity = "northwest") 

p_map_sampling <-
  c(img1_labeled, img2_labeled) %>% 
  image_append() 

plot(p_map_sampling)

image_write(p_map_sampling, path = "fig/map_sampling.png", format = "png", density = 600, flatten = TRUE)
