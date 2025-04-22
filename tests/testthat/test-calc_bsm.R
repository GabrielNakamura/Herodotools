library(Herodotools)

data("akodon_sites")
data("akodon_newick")


site_xy <- akodon_sites %>% 
  dplyr::select(LONG, LAT) 

akodon_pa <- akodon_sites %>% 
  dplyr::select(-LONG, -LAT)

spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
akodon_pa_tree <- akodon_pa[, spp_in_tree]


regions <-  
  Herodotools::calc_evoregions(
    comm = akodon_pa_tree,
    phy = akodon_newick
  )

site_region <- regions$cluster_evoregions

a_region <- Herodotools::get_region_occ(comm = akodon_pa_tree, site.region = site_region)


# ancestral reconstruction
load(file = system.file("extdata", "resDEC_akodon.RData", package = "Herodotools")) 


evoregion_df <- data.frame(
  site_xy, 
  site_region
)

# converting numbers to character
biogeo_area <- data.frame(biogeo = chartr("12345", "ABCDE", evoregion_df$site_region)) 

# getting the ancestral range area for each node 
node_area <- 
  Herodotools::get_node_range_BioGeoBEARS(
    resDEC,
    phyllip.file = here::here("inst", "extdata", "geo_area_akodon.data"),
    akodon_newick,
    max.range.size = 4 
  )


bsm_akodon <- calc_bsm(
  BioGeoBEARS.data = resDEC,
  phyllip.file = here::here("inst", "extdata", "geo_area_akodon.data"),
  tree.path = here::here("inst", "extdata", "akodon.new"),
  max.maps = 100, 
  n.maps.goal = 50,
  seed = 1234,
)


bsm_akodon$RES_clado_events_tables[[1]]



