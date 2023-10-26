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
    comm_names <- rownames(comm)
    list_res <- vector(mode = "list", length = nrow(comm))
    list_res2 <- vector(mode = "list", length = nrow(comm))
    names(list_res) <- rownames(comm)
    
    for(i in 1:length(list_res)){
      # processing results from reconstruction
      reconstruction <- phyloregion::dense2long(t(ifelse(ada.obj$reconstruction >= threshold, 1, 0))) # nodes predicted from reconstruction
      community <- phyloregion::dense2long(ada.obj$phylogeny) # Nodes extracted from community phylogeny 
      nodes_reconstruction <-
        reconstruction %>% 
        subset(grids == comm_names[i]) # reconstruction
      nodes_comm <- 
        community %>% 
        subset(grids == comm_names[i]) # phylogeny
      
      # finding the most basal node among all set of nodes including reconstruction and community nodes for a given community
      node_sequence <- gsub(pattern = "N", replacement = "", unique(c(nodes_reconstruction$nodes, nodes_comm$nodes)))
      data_nodes <- data.frame(nodes = unique(c(nodes_reconstruction$nodes, nodes_comm$nodes)), node_sequence = as.numeric(node_sequence))
      basal_pos <- which(data_nodes$node_sequence == min(data_nodes$node_sequence))
      tree_potential <- ape::extract.clade(phy = tree_test, node = data_nodes[basal_pos, 1])
      
      # naming node categories
      # ISDiv - In Situ Diversification nodes
      nodes_insitu <- intersect(nodes_reconstruction$nodes, nodes_comm$nodes)
      tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)], "IS", sep = "_") # nomeando nos com IS
      
      # IM - Immigration nodes
      node_immigration <- setdiff(nodes_comm$nodes, nodes_reconstruction$nodes)
      tree_potential$node.label[match(node_immigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_immigration, tree_potential$node.label)], "IM", sep = "_")
      
      # EM - Emmigration nodes
      node_emmigration <- setdiff(nodes_reconstruction$nodes, nodes_comm$nodes)
      tree_potential$node.label[match(node_emmigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_emmigration, tree_potential$node.label)], "EM", sep = "_")
      
      # ESDiv - Ex Situ Diversification nodes
      node_exsitu <- tree_potential$node.label[-grep(pattern = "_", tree_potential$node.label)] 
      tree_potential$node.label[match(node_exsitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_exsitu, tree_potential$node.label)], "ESD", sep = "_")
      
      # Organizing data to calculate PD components
      table_tree_potential2 <- 
        table_tree_potential %>% 
        dplyr::mutate(pres = ifelse(label %in% colnames(comm_obs), "pres", "abs")) %>% 
       dplyr::mutate(ancestor = table_tree_potential$label[table_tree_potential$parent]) %>% 
       dplyr::mutate(descendant = table_tree_potential$label[table_tree_potential$node]) %>% 
       dplyr::mutate(ancestor1 = gsub(pattern = ".*_", replacement = "", x = ancestor)) %>% 
       dplyr::mutate(descendent1 = gsub(pattern = ".*_", replacement = "", x = descendant)) %>% 
       dplyr::mutate(partition.IS = ifelse(ancestor1 == "IS" & descendent1 == "IS" |
                                       ancestor1 == "IS" & descendent1 == "EM" |
                                       ancestor1 == "EM" & descendent1 == "IS" | 
                                       ancestor1 == "EM" & descendent1 == "EM" |
                                       ancestor1 == "IS" & pres == "pres" |
                                       ancestor1 == "EM" & pres == "pres",
                                     "IS", NA)) %>% 
        dplyr::mutate(partition.IM = ifelse(ancestor1 == "IM" & descendent1 == "IS" |
                                       ancestor1 == "IM" & descendent1 == "EM",
                                     "IM", NA)) %>% 
        dplyr::mutate(partition.EM = ifelse(ancestor1 == "IS" & descendent1 == "IM" |
                                       ancestor1 == "EM" & descendent1 == "ESD" |
                                       ancestor1 == "EM" & descendent1 == "IM" | 
                                       ancestor1 == "EM" & pres == "abs" & descendent1 != "IS", 
                                     "EM", NA)) %>% 
        dplyr::mutate(partition.ESD = ifelse(ancestor1 == "ESD" & descendent1 == "ESD" |
                                        ancestor1 == "ESD" & pres == "abs", "ESD", NA)) %>% 
        dplyr::mutate(partition.unknown = ifelse(ancestor1 == "IM" & descendent1 == "IM", "unknown", NA))
      
      # calculating PD components
      PDinsitu <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.IS == "IS") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = T)
      PDimmigration <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.IM == "IM") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = T)
      
      PDemmigration <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.EM == "EM") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = TRUE)
      PDexsitu <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.ESD == "ESD") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = TRUE)
      
      PDunknown <- 
        table_tree_potential2 %>% 
        filter(partition.unknown == "unknown") %>% 
        select(branch.length) %>% 
        sum(na.rm = TRUE)
      
      
      PDtotal <- PDinsitu + PDimmigration + PDemmigration + PDexsitu + PDunknown
      
      #joining all results
      data_res <- 
        data.frame(partition = c("PDinsitu", "PDimmigration", "PDemmigration", "PDexsitu", "PDtotal", "PDunknown"), 
                   value = c(PDinsitu, PDimmigration, PDemmigration, PDexsitu, PDtotal, PDunknown), community = comm_names[i])
      table_tree_potential_res <- table_tree_potential2
      list_res[[i]] <- data_res
      list_res2[[i]] <- table_tree_potential_res
    }
    df_res <- do.call(rbind, list_res)
    list_res_final <- vector(mode = "list")
    list_res_final$PD_decomposition <- df_res
    list_res_final$tree_table_potential <- list_res2
    return(list_res_final)
  }