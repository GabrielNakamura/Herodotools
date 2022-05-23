
# reading data and libraries ----------------------------------------------

library(BioGeoBEARS)
library(devtools)
library(dplyr)
library(magrittr)
load_all() # temporary - load package functions

# Data

ancestral_rec <- readRDS(file = here("examples", "output", "resDEC_akodon.rds")) # ancestral reconstruction
akodon.sites <- read.table(
  here("inst", "extdata", "Table_Akodon_coords_pa.txt"),
  header = T
)
akodon.newick <- read.tree(here("examples", "data", "akodon.new"))

## Processing data

site.xy <- akodon.sites %>% 
  dplyr::select(LONG, LAT)

akodon.pa <- akodon.sites %>% 
  dplyr::select(-LONG, -LAT)

spp.in.tree <- names(akodon.pa) %in% akodon.newick$tip.label

akodon.pa.tree <- akodon.pa[, spp.in.tree]


# getting ancestral node data  ------------------------------------------------------

node.area <- 
  get.node.range_BioGeoBEARS(
    ancestral_rec,
    phyllip.file = here("examples", "data", "geo_area_akodon.data"),
    akodon.tree,
    max.range.size = 4 
  )
node_data <- data.frame(area = node.area[, 2])
nodes.biomes <- data.frame(area_node = node.area$biome[-c(1:30)])
rownames(nodes.biomes) <- paste("N", node.area[-c(1:30), 1], sep = "")



# calculating age arrival -------------------------------------------------

biogeo_area <- data.frame(biogeo = chartr("12345", "ABCDE", evoregion.df$site.region)) # converting numbers to character
age_comm <- age_arrival(W = akodon.pa, tree = akodon.pa.tree, ancestral.area = nodes.biomes, biogeo = biogeo_area) # calculating age arrival 
?age_arrival
