#' Phylogenetic diversity decomposition 
#'
#' @param comm a community matrix with sites in rows and species in columns
#' @param ada.obj an object calculated from \code{\link{ada_core}} containing community reconstruction 
#' @param phy a newick object containing phylogenetic relationships among species
#' @param threshold a scalar indicating the threshold used to consider the presence of a species in a community
#'
#' @return a list of length two. The element of the list named decomposition is a data frame with decomposition values of PD for all communities.
#'     The element of the list named tree_table_potential is a tibble with node information for the potential tree for each community
#' @export
#'
#' @examples
PD_decomposition <- 
  function(comm, ada.obj, phy, threshold){
    comm <- ifelse(comm >= 1, 1, 0)
    phy <- ape::makeNodeLabel(phy = phy, method = "number", prefix = "Node")
    comm_names <- rownames(comm)
    list_res <- vector(mode = "list", length = nrow(comm)) # object to receive potential tree partition
    list_res2 <- vector(mode = "list", length = nrow(comm)) # object to receive potential tree partition
    list_res3 <- vector(mode = "list", length = nrow(comm)) # object to receive Faith partition
    list_res4 <- vector(mode = "list", length = nrow(comm)) # object to receive Faith partition
    names(list_res) <- rownames(comm)
    
    reconstruction <- phyloregion::dense2long(t(ifelse(ada.obj$reconstruction >= threshold, 1, 0))) # nodes predicted from reconstruction
    community <- phyloregion::dense2long(ada.obj$phylogeny) # Nodes extracted from community phylogeny 
    phy_tibble <- tidytree::as_tibble(phy)
    pb <- txtProgressBar(min = 0,      
                         max = length(list_res), 
                         style = 3,    
                         width = 50,   
                         char = "=")   
    for(i in 1:length(list_res)){
      # processing results from reconstruction
      nodes_reconstruction <-
        reconstruction %>% 
        subset(grids == comm_names[i]) # reconstruction
      
      nodes_comm <- 
        community %>% 
        subset(grids == comm_names[i]) # phylogeny
      spp_comm <- names(which(comm[i, ] == 1)) # species in community
      comm_tible <- phy_tibble[match(c(spp_comm, nodes_comm$species), phy_tibble$label), ]
      offspring_data <- tidytree::offspring(phy_tibble, nodes_reconstruction$species)
      if(class(offspring_data)[1] == "list"){
        df_offspring_rec <- do.call(rbind, tidytree::offspring(phy_tibble, nodes_reconstruction$species))
      } else{
        df_offspring_rec <- offspring_data
      }
      
      
      # starting calculation of PD components
      if(is.null(df_offspring_rec) == TRUE){ # when communities were not reconstructed at all
        tree_potential <- ape::keep.tip(phy = phy, tip = spp_comm)
      } else{
        spp_potential <- 
          df_offspring_rec %>% 
          distinct(label) # species in observed communities
        spp_potential2 <- phy$tip.label[phy$tip.label %in% spp_potential$label] # species names from nodes estimated in reconstruction
        spp_potential_all <- unique(c(spp_potential2, spp_comm)) # joining community and reconstruction 
        tree_potential <- ape::keep.tip(phy = phy, tip = spp_potential_all) # keeping only species observed in communities and estimated in reconstruction
      }
      
      # naming node categories
      # ISDiv - In Situ Diversification nodes
      nodes_insitu <- intersect(nodes_reconstruction$species, nodes_comm$species)
      tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)], "IS", sep = "_") # nomeando nos com IS
      
      # IM - Immigration nodes
      node_immigration <- setdiff(nodes_comm$species, nodes_reconstruction$species)
      tree_potential$node.label[match(node_immigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_immigration, tree_potential$node.label)], "IM", sep = "_")
      
      # EM - emigration nodes
      node_emigration <- setdiff(nodes_reconstruction$species, nodes_comm$species)
      tree_potential$node.label[match(node_emigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_emigration, tree_potential$node.label)], "EM", sep = "_")
      
      # ESDiv - Ex Situ Diversification nodes
      node_exsitu <- tree_potential$node.label[-grep(pattern = "_", tree_potential$node.label)] 
      tree_potential$node.label[match(node_exsitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_exsitu, tree_potential$node.label)], "ESD", sep = "_")
      
      table_tree_potential <- tidytree::as_tibble(tree_potential)
      tree_faith <- ape::keep.tip(phy = tree_potential, tip = spp_comm) # keeping only species observed in communities and estimated in reconstruction
      table_tree_faith <- tidytree::as_tibble(tree_faith)
      
      comm_faith <- matrix(rep(1, length(spp_comm)), nrow = 1, ncol = length(spp_comm), 
                           dimnames = list(comm_names[i], spp_comm))
      comm_potential <- matrix(rep(1, length(spp_potential_all)), nrow = 1, ncol = length(spp_potential_all),
                               dimnames = list(comm_names[i], spp_potential_all))
      
      # Organizing data to calculate PD components
      partition_potential <-  
        suppressMessages(
          suppressWarnings(
            partition_pd(data = table_tree_potential, comm = comm_potential, phy = phy)
          )
        ) 
      partition_faith <- 
        suppressMessages(
          suppressWarnings(
            partition_pd(data = table_tree_faith, comm = comm_faith, phy = phy)
          )
        )
      
      # results partition potential
      list_res[[i]] <- partition_potential$data_res
      list_res2[[i]] <- partition_potential$table_tree 
      
      # results partition faith tree
      list_res3[[i]] <- partition_faith$data_res
      list_res4[[i]] <- partition_faith$table_tree
      
      # progress bar
      setTxtProgressBar(pb, i)
    }
    
    # summarizing results
    decomp_potential <- organize_pd_decomp(list1.decomp = list_res, list2.treetable = list_res2) # function used to organize the results from pd decomposition
    decomp_faith <- organize_pd_decomp(list1.decomp = list_res3, list2.treetable = list_res4)
    list_res_final <- vector(mode = "list")
    list_res_final$decomp_potential <- decomp_potential
    list_res_final$decomp_faith <- decomp_faith
    return(list_res_final)
  }