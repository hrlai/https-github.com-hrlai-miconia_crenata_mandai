
library(tidyverse)
library(brms)
library(sf)



# Data --------------------------------------------------------------------

# remove transect C63 for now because I couldn't find its spatial data
quadrat  <- 
  read_csv("data/quad.csv") %>% 
  filter(Tran != "C63")
transect <- 
  read_csv("data/spatialsurveyed.csv")%>% 
  filter(Tran != "C63")

# canopy opening
canopy <- 
  read_csv("data/canopy.csv") %>% 
  select(Tran = `Transect ID`,
         Quad = `Quadrat ID (Distance)`,
         Quad_branch = `Quadrat branching`,
         Densiometer = `Densiometer reading`,
         Densio_dir = Direction) %>% 
  mutate(Densiometer = as.numeric(Densiometer) / 96,
         Densiometer = (Densiometer * (n()-1) + 1/2) / n()) %>% 
  filter(Tran != "C63",
         !is.na(Densiometer)) %>% 
  arrange(Tran) %>% 
  # summarise mean and sd to fed into measurement error model
  group_by(Tran, Quad, Quad_branch) %>% 
  summarise(Densio    = mean(Densiometer, na.rm = TRUE),
            Densio_sd = sd(Densiometer))

# leaf litter depth
litter <- 
  quadrat %>% 
  select(Tran, Quad, Litter) %>% 
  left_join(transect) %>% 
  arrange(Tran) %>% 
  filter(!is.na(Litter)) 

# region of interest
roi <- 
  # st_read("data/GIS/ccnr/ccnr.shp") %>% 
  st_read("data/Survey area.shp") %>% 
  st_transform(crs = st_crs("EPSG:32648")) 
bbox_crop <- st_bbox(roi)
bbox_crop["ymin"] <- 153000
roi <- st_crop(roi, bbox_crop)

# water edges (reservoir)
reservoir <- 
  st_read("data/Reservoir.shp") %>% 
  st_transform(crs = st_crs("EPSG:32648"))

# roads 
roads <- 
  st_read("data/Updated roads.shp") %>% 
  st_transform(crs = st_crs("EPSG:32648"))

# sampled coords
coord <- 
  quadrat %>% 
  select(Patch, Tran, Quad, Qb, CorrN, CorrE) %>% 
  st_as_sf(coords = c("CorrE", "CorrN"),
           crs = st_crs("EPSG:32648")) %>% 
  st_geometry()

# distance to nearest road
dist_road <- 
  coord %>% 
  st_distance(st_combine(roads))

# distance to nearest water edge
dist_reservoir <- 
  coord %>% 
  st_distance(st_combine(reservoir))

# distance to nearest vegetation edge
veg_boundary <- 
  st_read("data/veg_edge_edit.shp") %>% 
  st_transform(crs = st_crs("EPSG:32648"))
dist_veg <- 
  coord %>% 
  st_distance(st_combine(veg_boundary))

# distance to edge (can be either vegetation or reservoir)
dist_edge <- 
  coord %>% 
  st_distance(st_cast(roi, "MULTILINESTRING") %>% st_combine)

# calculate road density again
# Hao Ran setting the radius to 250m after eyeballing the distribution of 
# distance from edge
# this is slow, only rerun when necessary
rerun <- FALSE
if (rerun) {
  rd_250_buffer <- 
    coord %>% 
    st_buffer(250) %>% 
    # clip against reservoir to remove crossing water bodies
    st_difference(reservoir) %>% 
    st_geometry()
  # now we clip a (non-contiguous) buffer if it crosses to the opposite shore, so we don't 
  # measure roads on the other side
  rd_250_buffer_clipped <- rd_250_buffer
  for (i in seq_along(rd_250_buffer_clipped)) {
    buffers_tmp <- st_cast(rd_250_buffer[i], "POLYGON")
    buffer_select <- as.numeric(st_within(coord[i], buffers_tmp))
    rd_250_buffer_clipped[i] <- buffers_tmp[buffer_select]
  }
  # annoyingly, st_intersection doesn't return zero when there is no intersection
  # so we need some manual looping
  rd_250 <- numeric(length(rd_250_buffer_clipped))
  for (i in 1:length(rd_250_buffer_clipped)) {
    out <- 
      st_intersection(x = st_combine(roads), 
                      y = rd_250_buffer_clipped[i]) %>% 
      st_length() %>% 
      as.numeric()
    if (length(out) > 0)  { rd_250[i] <- out }
  }
  write_rds(rd_250, "data/rd_250.rds")
}
rd_250 <- read_rds("data/rd_250.rds")

# Clidemia occupancy
clidemia <-
  quadrat %>% 
  # join road density at 250 m radius
  mutate(RD250 = as.numeric(rd_250 / (pi*(250^2)))) %>% 
  # join new distance measures
  mutate(dist_veg = as.numeric(dist_veg), 
         dist_road = as.numeric(dist_road),
         dist_edge = as.numeric(dist_edge),
         dist_reservoir = as.numeric(dist_reservoir)) %>% 
  filter(!is.na(Clid),
         !is.na(Litter)) %>% 
  # join densiometer measurement
  left_join(canopy) %>%
  select(Patch, Tran, 
         Northing = CorrN, Easting = CorrE, Dir,
         Quad, Qb, 
         Clid,
         Litter,
         Densio, Densio_sd,
         edgetype,
         distedge, dist_edge,
         dist_veg, dist_road, dist_reservoir,
         RD200, RD250, RD400, RD800) %>% 
  arrange(Tran)

# scale predictors
clidemia_s <- 
  clidemia %>% 
  mutate_at(vars(distedge, dist_edge,
                 dist_veg, dist_road, dist_reservoir,
                 RD200, RD250, RD400, RD800),
            ~ as.numeric(scale(.))) %>%
  mutate_at(vars(Northing, Easting),
            ~ . / 1000) %>% 
  mutate(Litter = Litter / 100)




# Model -------------------------------------------------------------------

bform_veg <-
  bf(Clid ~ 1 + Litter + mi(Densio) + dist_veg + s(Easting, Northing, k = 50),
     center = TRUE) +
  bf(Densio | mi() ~ 1 + dist_veg + s(Easting, Northing, k = 50)) +
  bf(Litter ~ 1 + dist_veg + s(Easting, Northing, k = 50)) +
  set_rescor(FALSE)

fit_veg <-
  brm(bform_veg,
      data = clidemia_s,
      family =
        list(bernoulli(),
             Beta(),
             hurdle_lognormal()),
      file = "out/mod_brms_dist_veg.rds",
      cores = 4)
