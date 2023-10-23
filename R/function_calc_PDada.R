
#' Calculation of PD using ada object for one community
#'
#' @param ada.obj 
#' @param threshold 
#' @param comm.name 
#'
#' @return
#' @export
#'
#' @examples
calc_PD_ada <- 
  function(ada.obj, threshold = 0.6, comm.name){
    reconstruction <- phyloregion::dense2long(t(ifelse(ada.obj$reconstruction >= 0.6, 1, 0))) # nodes predicted from reconstruction
    phylo_only <- phyloregion::dense2long(ada.obj$phylogeny) # Nodes extracted from phylogeny 
    
    # calculation for each community
    nodes_reconstruction <-
      reconstruction %>% 
      subset(grids == comm.name) # reconstruction
    
    nodes_phylogeny <- 
      phylo_only %>% 
      subset(grids == comm.name) # phylogeny
    
    # nodes community
    
    comm_obs_total <- which(x["comm_10", ] == 1) # spp observed
    
    # core components
    phy <- ape::makeNodeLabel(phy = tree, method = "number", prefix = "Node")
    phy_base <- as(phy, "phylo4")
    # tip_total <- phylobase::descendants(phy = phy_base, node = unique(c(nos_reconstrucao$species, nos_comm$species))) # phylogeny and reconstruction
    # names_total <- unique(unlist(lapply(tip_total, names))) # species names phylogeny + reconstruction
    nodes_potential <- unique(c(nodes_reconstruction$species, nodes_phylogeny$species))
    
    # observed PD
    node_local <- intersect(nodes_reconstruction$species, nodes_phylogeny$species) # nodes with in situ diversification
    tip_local <- phylobase::descendants(phy = phy_base, node = node_local) # tips from local nodes
    names_local <- unique(unlist(lapply(tip_local, names))) # potential in situ diversification 
    tip_obs <- names(comm_obs_total)
    tips_in_situ <- intersect(names_local, tip_obs) # in situ species
    tip_immigration <- setdiff(tip_obs, tips_in_situ)
    tip_potential <- phylobase::descendants(phy = phy_base, node = nodes_potential)
    names_potential <- unique(unlist(lapply(tip_potential, names)))
    
    # PDs
    PD_obs_base <- picante::pd(samp = t(as.matrix(x[comm.name, ])), tree = phy, include.root = TRUE)
    PD_insitu <- 
      picante::pd(samp = matrix(1, nrow = 1, ncol = length(tips_in_situ), dimnames = list("comm", tips_in_situ)),
                  tree = phy)
    PD_immigration <- picante::pd(samp = matrix(1, nrow = 1, ncol = length(tip_immigration), dimnames = list("comm", tip_immigration)),
                                  tree = phy)
    PD_obs <- PD_insitu[1, 1] + PD_immigration[1, 1]
    
    PDpotencial <-  picante::pd(samp = matrix(1, nrow = 1, ncol = length(names_potential), dimnames = list("comm", names_potential)),
                                tree = phy)
    
    res_PD <- rbind(obs = PD_obs_base, in.situ = PD_insitu, immigration = PD_immigration, potential = PDpotencial)
    res_PD$comm <- comm.name
    return(res_PD)
  }

