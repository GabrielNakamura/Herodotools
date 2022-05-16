
# load packages -----------------------------------------------------------

library(devtools)
library(tidyverse)
library(ape)
library(here)
library(rnaturalearth)
library(viridis)
library(devtools)
library(terra)
library(sf)
library(patchwork)
load_all()

# load data ---------------------------------------------------------------

akodon.tree <- read.nexus(
   here("inst", "extdata", "tree_akodon.nexus")  
)

akodon.sites <- read.table(
  here("inst", "extdata", "Table_Akodon_coords_pa.txt"),
  header = T
)

site.xy <- akodon.sites %>% 
  dplyr::select(LONG, LAT)

akodon.pa <- akodon.sites %>% 
  dplyr::select(-LONG, -LAT)

# |- exploring data ----

# |-- names macthing ----

akodon.tree$tip.label <- akodon.tree$tip.label %>% 
  str_replace_all("Akodon", "A")

### save a newick to use in BioGeoBEARS
write.tree(akodon.tree, here("examples", "data", "akodon.new"))

akodon.newick <- read.tree(here("examples", "data", "akodon.new"))

spp.in.tree <- names(akodon.pa) %in% akodon.newick$tip.label

akodon.pa.tree <- akodon.pa[, spp.in.tree]


# |-- spatial pattern ----

coastline <- ne_coastline(returnclass = "sf")
map.limits <- list(
  x = c(-95, -30),
  y = c(-55, 12)
)

ggplot(site.xy) + 
  geom_raster(aes(x = LONG, y = LAT)) + 
  geom_sf(data = coastline) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme_bw()


# |-- richness ----

rich <- rowSums(akodon.pa)
rich.tree <- rowSums(akodon.pa.tree)

(map_rich <- 
bind_cols(site.xy, rich =  rich, rich.tree = rich.tree) %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = rich)) + 
  scale_fill_viridis(option = "plasma") +
  geom_sf(data = coastline) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme_bw()
)

(map_rich_tree <- 
bind_cols(site.xy, rich =  rich, rich.tree = rich.tree) %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = rich.tree)) + 
  scale_fill_viridis(option = "plasma") +
  geom_sf(data = coastline) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme_bw()
)

# evoregions --------------------------------------------------------------

set.seed(12)
regions <- evoregions(
  comm = akodon.pa.tree,
  phy = akodon.tree, 
  max.n.clust = 10)

site.region <- regions$Cluster_Evoregions

evoregion.df <- data.frame(
  site.xy, 
  site.region
)

r.evoregion <- rast(evoregion.df)

sf.evoregion <- as.polygons(r.evoregion) %>% 
  st_as_sf()

st_crs(sf.evoregion) <- st_crs(coastline)

col_five_hues <- c(
  "#3d291a",
  "#a9344f",
  "#578a5b",
  "#83a6c4",
  "#fcc573"
)

# |- map evoregion ----
(map_evoregion <- 
evoregion.df %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = site.region)) + 
  scale_fill_manual(
    name = "Evoregions", 
    labels = LETTERS[1:5],
    values = rev(col_five_hues)
    ) +
  geom_sf(data = coastline) +
  geom_sf(
    data = sf.evoregion, 
    color = "#040400",
    fill = NA, 
    size = 0.2) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  ggtitle("Evoregion for Akodon Genus") + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    plot.title.position =  "plot"
  )
)

ggsave(
  here("examples", "output", "fig", "fig_01_akodon_evoregion.png"),
  map_evoregion,
  width = 4,
  height = 5
)

# |- evoregion's transition zones ----

# only axis with more than 5% of explained variance
axis_sel <- which(regions$PCPS$prop_explainded >= regions$PCPS$tresh_dist)
PCPS_thresh <- regions$PCPS$vectors[, axis_sel] 

# distance matrix from 4 significant PCPS axis
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# affiliation values for each assemblage 
afi <- afilliation.evoreg(phylo.comp.dist = dist_phylo_PCPS,
                          groups = regions$Cluster_Evoregions) 

sites <- bind_cols(site.xy, site.region =  site.region, afi)

# |- mapping evoregions and afilliation --------------------------------------

(map_afiliation <- 
sites %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = afilliation)) + 
  scale_fill_gradient(
    name = "Afiliation",
    low = "#8EC1CD",
    high = "#3B3C68") + 
  geom_sf(data = coastline) +
  geom_sf(
    data = sf.evoregion, 
    color = "#040400",
    fill = NA, 
    size = 0.2) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank()
  )
)


map_evoregion_afilliation <- map_evoregion + map_afiliation

ggsave(
  here("examples", "output", "fig", "fig_02_akodon_evo_affiliation.png"),
  map_evoregion_afilliation,
  width = 8,
  height = 5
)

# defining regional merbership of species ----------------------------------

akodon.evoregion.data <- 
  bind_cols(
    akodon.pa.tree,
    site.region =  site.region
  ) %>% 
  pivot_longer(
    cols = 1:30, 
    names_to = "species",
    values_to = "presence"
  ) %>% 
  filter(presence == 1) %>% 
  group_by(species, site.region) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(
    species.total = sum(n), 
    prec.occupation = n/species.total
  ) %>% 
  ungroup() %>% 
  filter(
    prec.occupation >= .25
  ) %>% 
  mutate(
    area = LETTERS[site.region]
  ) %>% 
  dplyr::select(species, area) %>% 
  mutate( value = 1) %>% 
  pivot_wider(
    id_cols = species,
    names_from = area, 
    names_sort = T,
    values_from = value,
    values_fill = 0
  ) %>% 
  as.data.frame()
  
species.names <- akodon.evoregion.data[,1]
a.regions <- akodon.evoregion.data[,-1]

rownames(a.regions) <- species.names

# save phyllip file
tipranges_to_BioGeoBEARS(
  a.regions, 
  filename = here("examples", "data", "geo_area_akodon.data"),
  areanames = NULL
  )


# run DEC model with BioGeoBears ------------------------------------------

## this run DEC model for akodon species using BioGeoBEARS
source(here("examples", "script", "e_01_run_DEC_model.R"))

## the object with results is 'resDEC'
saveRDS(resDEC, here("examples", "output", "resDEC_akodon.rds"))

# |- exploring DEC results ----
resDEC <- readRDS(here("examples", "output", "resDEC_akodon.rds"))

node.area <- 
get.node.range_BioGeoBEARS(
  resDEC,
  phyllip.file = here("examples", "data", "geo_area_akodon.data"),
  akodon.tree,
  max.range.size = 4 
)

nodes.biomes <- node.area$biome[-c(1:30)]

tip.biomes <- apply(a.regions, 1, function(x){
  paste(names(x)[x == 1], collapse = "")
})

order.tip <- match(
  akodon.newick$tip.label,
  rownames(a.regions)
  )


plot(akodon.newick)
nodelabels(nodes.biomes, cex = 0.75)
tiplabels(tip.biomes[(order.tip)], cex = 0.75)
axisPhylo()






