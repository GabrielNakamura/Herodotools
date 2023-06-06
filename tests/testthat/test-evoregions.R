## data for examples ----

data(akodon_sites) # occurrence data 
data(akodon_newick) # phylogenetic tree

akodon_pa <- akodon_sites %>% 
  dplyr::filter(
    dplyr::between(LAT, -30, -10),
    dplyr::between(LONG, -80, -60)
    ) %>%
  dplyr::select(-LONG, -LAT)

akodon_pa <- akodon_pa[,colSums(akodon_pa) != 0]

site_xy <- akodon_sites %>% 
  dplyr::filter(
    dplyr::between(LAT, -30, -10),
    dplyr::between(LONG, -80, -60)
    ) %>%
  dplyr::select(LONG, LAT)

phy_comm <- picante::match.phylo.comm(akodon_newick, akodon_pa)

phy <- phy_comm$phy
comm <- phy_comm$comm



test_that("seed for the names of the groups", {
  seed = 265
  max.n.clust = 6
  
  regions <- calc_evoregions(
    comm = comm, 
    phy = phy, 
    max.n.clust = max.n.clust,
    seed = seed)
  
  site_region1 <- regions$cluster_evoregions 
  
  regions <- calc_evoregions(
    comm = comm, 
    phy = phy, 
    max.n.clust = max.n.clust,
    seed = seed)
  
  site_region2 <- regions$cluster_evoregions 
  
  expect_equal(site_region1, site_region2)
})


