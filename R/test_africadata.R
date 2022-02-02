devtools::install_github("GabrielNakamura/Rrodotus", force = TRUE)

library(Rrodotus)
library(phyloregion)
data(africa)
sparse_comm <- africa$comm

# pb-sim as distance matrix -----------------------------------------------

tree <- africa$phylo
tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
pb <- phylobeta(sparse_comm, tree)
optimal <- optimal_phyloregion(x = pb$phylo.beta.sim, method = "average", k = 20)
phylo_regionalization <- phyloregion(x = pb$phylo.beta.sim, k = optimal$optimal$k, method = "average")

# phylogenetic fuzzy as distance matrix -----------------------------------

comm_dense <- phyloregion::sparse2dense(sparse_comm)
matrix_p_nigeria <- SYNCSA::matrix.p(comm = comm_dense, phylodist = cophenetic(tree))
beta_fuzzy_bray <- sqrt(vegan::vegdist(matrix_p_nigeria$matrix.P, method = "bray"))
optimal_fuzzy <- phyloregion::optimal_phyloregion(x = beta_fuzzy_bray, method = "average", k = 20)
phylo_regionalization_fuzzy <- phyloregion::phyloregion(x = beta_fuzzy_bray, 
                                                        k = optimal_fuzzy$optimal$k,
                                                        method = "average")



# classification with evoregion -------------------------------------------

phylo_evoregion <- Rrodotus::evoregions(comm = comm_dense, phy = phytools::force.ultrametric(tree), max.n.clust = NULL,
           max.n.clust.method = "elbow",
           method.dist = "bray",
           tresh.dist = 0.05, 
           method.clust = "kmeans",
           stat.clust = "BIC", 
           n.iter.clust = 1e7, 
           criterion.clust = "diffNgroup")

grp_evoregions <- phylo_evoregion$Cluster_Evoregions$grp
df_evoregion <- data.frame(grids = names(phylo_evoregion$Cluster_Evoregions$grp), 
                           cluster = phylo_evoregion$Cluster_Evoregions$grp)
sf_evoregion <- sf_africa %>%
  st_transform(crs = "+proj=robin") %>% 
  left_join(df_evoregion)


# calculating affiliation of each cell ------------------------------------




help(package = "Rrodotus")
afilliation_evoregion <- afilliation.evoreg(evo.vectors = phylo_evoregion, method = "euclidean")
afilliation_belonging <- affiliation_evoregion[[2]]
affiliation_belonging <- data.frame(afilliation = afilliation_belonging[, 1], 
                                    group = afilliation_belonging[, 2],
                                    grids = rownames(afilliation_belonging))
sf_evoregion_belonging <- sf_evoregion %>%
  left_join(affiliation_belonging)

ggplot() +
  geom_sf(data = sf_africa, aes(geometry = geometry), fill = NA) +
  geom_sf(data = sf_evoregion_belonging %>% filter(group == 4), aes(geometry = geometry, 
                                     fill = afilliation),
          size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
    geom_sf(data = sf_evoregion_belonging %>% filter(group == 1), aes(geometry = geometry, 
                                                                    fill = afilliation),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "TealGrn", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
  geom_sf(data = sf_evoregion_belonging %>% filter(group == 5), aes(geometry = geometry, 
                                                                    fill = afilliation),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "Magenta", 
                                  direction = 1, 
                                  limits = c(0, 1)) 


quartz()
rcartocolor::display_carto_all(
)


# spatialization ----------------------------------------------------------

library(ggplot2)  
library(tidyverse)
library(sf)
library(patchwork)

sf_africa <- sf::st_as_sf(africa$polys)
sf_phyloregion <- sf_africa %>%
  st_transform(crs = "+proj=robin") %>% 
  left_join(phylo_regionalization$region.df)

map_phyloregion <- ggplot() +
  geom_sf(data = sf_phyloregion, aes(geometry = geometry, 
                                   fill = as.factor(cluster)),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_d(type = "qualitative",
                                  palette = 6) +
  guides(fill = guide_colorbar(barheight = unit(2.3, units = "mm"),  
                               barwidth = unit(100, units = "mm"),
                               direction = "horizontal",
                               ticks.colour = "grey20",
                               title.position = "top",
                               label.position = "bottom",
                               title.hjust = 0.5)) 


sf_fuzzy <- sf_africa %>%
  st_transform(crs = "+proj=robin") %>% 
  left_join(phylo_regionalization_fuzzy$region.df)

map_phyloregion_fuzzy <- ggplot() +
  geom_sf(data = sf_fuzzy, aes(geometry = geometry, 
                                     fill = as.factor(cluster)),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_d(type = "qualitative",
                                  palette = 6) +
  guides(fill = guide_colorbar(barheight = unit(2.3, units = "mm"),  
                               barwidth = unit(100, units = "mm"),
                               direction = "horizontal",
                               )
  )

map_evoregion <- ggplot() +
  geom_sf(data = sf_evoregion, aes(geometry = geometry, 
                               fill = as.factor(cluster)),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_d(type = "qualitative",
                                  palette = 6) +
  guides(fill = guide_colorbar(barheight = unit(2.3, units = "mm"),  
                               barwidth = unit(100, units = "mm"),
                               direction = "horizontal",
  )
  )

map_phyloregion + map_phyloregion_fuzzy + map_evoregion
