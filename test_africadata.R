# devtools::install_github("GabrielNakamura/Rrodotus", force = TRUE)
help(package = "Rrodotus")

# data and packages -------------------------------------------------------

# library
library(Rrodotus)
library(phyloregion)
library(magrittr)
library(ggplot2)  
library(tidyverse)
library(sf)
library(patchwork)


# data
data(africa)
sparse_comm <- africa$comm
sf_africa <- sf::st_as_sf(africa$polys)


# pb-sim as distance matrix -----------------------------------------------

tree <- africa$phylo
tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
pb <- phylobeta(sparse_comm, tree)
optimal <- optimal_phyloregion(x = pb$phylo.beta.sim, method = "average", k = 20)
phylo_regionalization <- phyloregion(x = pb$phylo.beta.sim,
                                     k = optimal$optimal$k,
                                     method = "average",
                                     shp = africa$polys)

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
  sf::st_transform(crs = "+proj=robin") %>% 
  dplyr::left_join(df_evoregion)


# calculating affiliation of each cell ------------------------------------

afilliation_evoregion <- afilliation.evoreg(evo.classification = phylo_evoregion, method = "euclidean")
afilliation_phyloregion <- afilliation.evoreg(evo.classification = phylo_regionalization, distance = pb$phylo.beta.sim)
NA_rm_phyloregion <- which(is.na(match(rownames(afilliation_phyloregion), rownames(afilliation_evoregion))) == TRUE) # check why NA in phyloregion
afilliation_phyloregion <- afilliation_phyloregion[-NA_rm_phyloregion, ]
afilliation_evoregion <- afilliation_evoregion[match(rownames(afilliation_phyloregion), rownames(afilliation_evoregion)), ]



# spatialization of regions and afilliation -------------------------------



afilliation_grid <- data.frame(afilliation_evoregion = afilliation_evoregion[, 1], 
                               group_evoregion = afilliation_evoregion[, 2],
                               afilliation_phyloregion = afilliation_phyloregion[, 1],
                               group_phyloregion = afilliation_phyloregion[, 2],
                               grids = rownames(afilliation_phyloregion))

sf_regionalization_membership <- sf_africa %>%
  dplyr::left_join(afilliation_grid) %>% 
  dplyr::left_join(phylo_regionalization$region.df[, c(1,4)])

ggplot() +
  geom_sf(data = sf_regionalization_membership,
          aes(geometry = geometry, 
              fill = as.factor(group_phyloregion)
              )
          )



ggplot() + 
  geom_sf(data = sf_regionalization_membership, aes(fill = as.factor(group_phyloregion), alpha = afilliation_phyloregion, colour = NA)) + 
  scale_fill_manual(values = levels(as.factor(sf_regionalization_membership$COLOURS)), alpha = sf_regionalization_membership$afilliation_phyloregion)

ggplot(data = sf_regionalization_membership) +
  geom_sf(aes(geometry = geometry, fill = afilliation_evoregion)) +
  scale_fill_manual(values = levels(as.character(sf_regionalization_membership$COLOURS)))


ggplot() +
  geom_sf(data = sf_regionalization_membership, aes(geometry = geometry), fill = NA, labels = NA) +
  geom_sf(data = sf_regionalization_membership %>% dplyr::filter(group_phyloregion == 1), aes(geometry = geometry, 
                                     fill = afilliation_phyloregion),
          size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "Peach", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
    geom_sf(data = sf_regionalization_membership %>% dplyr::filter(group_phyloregion == 2), aes(geometry = geometry, 
                                                                    fill = afilliation_phyloregion),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "TealGrn", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
  geom_sf(data = sf_regionalization_membership %>% dplyr::filter(group_phyloregion == 3), aes(geometry = geometry, 
                                                                    fill = afilliation_phyloregion),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "Magenta", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
  geom_sf(data = sf_regionalization_membership %>% dplyr::filter(group_phyloregion == 4), aes(geometry = geometry, 
                                                                                              fill = afilliation_phyloregion),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "Mint", 
                                  direction = 1, 
                                  limits = c(0, 1)) +
  ggnewscale::new_scale("fill") +
  geom_sf(data = sf_regionalization_membership %>% dplyr::filter(group_phyloregion == 5), aes(geometry = geometry, 
                                                                                              fill = afilliation_phyloregion),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "Teal", 
                                  direction = 1, 
                                  limits = c(0, 1)) 

quartz()

rcartocolor::display_carto_all()

phyloregion::plot_swatch


# spatialization ----------------------------------------------------------


sf_phyloregion <- sf_africa %>%
  st_transform(crs = "+proj=robin") %>% 
  dplyr::left_join(phylo_regionalization$region.df)

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
  dplyr::left_join(phylo_regionalization_fuzzy$region.df)

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
