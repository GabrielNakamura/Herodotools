
# load packages -----------------------------------------------------------

library(devtools)
library(tidyverse)
library(ape)
library(here)
library(rnaturalearth)
library(viridis)
library(devtools)
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
write.tree(akodon.tree, here("examples", "akodon.new"))

akodon.newick <- read.tree(here("examples", "akodon.new"))

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

regions <- evoregions(
  comm = akodon.pa.tree,
  phy = akodon.tree, 
  max.n.clust = 10)

length(regions$PCPS)
site.region <- regions$Cluster_Evoregions

col_five_hues <- c(
  "#3d291a",
  "#a9344f",
  "#578a5b",
  "#83a6c4",
  "#fcc573"
)

col_blue_gold_red <- c(
  "#303260",
  "#88beca",
  "#fcc573",
  "#7b5a28",
  "#EC7D2E",
  "#a9344f"
)

# |- map evoregion ----
(map_evoregion <- 
bind_cols(site.xy, site.region =  site.region) %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = site.region)) + 
  scale_fill_manual(
    name = "Evoregions", 
    labels = LETTERS[1:5],
    values = rev(col_five_hues)
    ) +
  geom_sf(data = coastline) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  ggtitle("Akodon's Evoregion") + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    plot.title.position =  "plot"
  )
)

ggsave(
  here("examples", "fig_01_akodon_evoregion.png"),
  map_evoregion,
  width = 4,
  height = 5
)

# |- evoregion's transition sozes ----
#### afiliation tem problemas!
afi <- afilliation.evoreg(regions)

sites <- bind_cols(site.xy, site.region =  site.region, afi)
sum(!sites$site.region == sites$group)
bind_cols(site.xy, site.region =  site.region, afi) %>% 
  ggplot() + 
  geom_raster(aes(x = LONG, y = LAT, fill = afilliation)) + 
  geom_sf(data = coastline) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


# definig regional merbership of species ----------------------------------

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
tipranges_to_BioGeoBEARS(a.regions, 
                         filename = here("examples", "geo_area_akodon.data"),
                         areanames = NULL)


# run DEC model with BioGeoBears ------------------------------------------

## this run DEC model for akodon species using BioGeoBEARS
source(here("examples", "e_01_run_DEC_model.R"))

## the object with results is 'resDEC'

# |- exploring DEC results ----
node.area <- 
get.node.range_BioGeoBEARS(
  resDEC,
  phyllip.file = "lagrange_area_data_file.data",
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



