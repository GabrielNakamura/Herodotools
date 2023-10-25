#' Phylogenetic diversity decomposition 
#'
#' @param comm a community matrix with sites in rows and species in columns
#' @param ada.obj an object calculated from \code{\link{ada_core}} containing community reconstruction 
#' @param phy a newick object containing phylogenetic relationships among species
#' @param threshold a scalar indicating the threshold used to consider the presence of a species in a community
#'
#' @return a data frame with three columns. partition indicates the PD component, value is the values of partitions, community is the name of community
#' @export
#'
#' @examples
PD_decomposition <- 
  function(comm, ada.obj, phy, threshold){
    comm <- ifelse(comm >= 1, 1, 0)
    comm_names <- rownames(comm)
    list_res <- vector(mode = "list", length = nrow(comm))
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
      table_tree_potential <- tidytree::as_tibble(tree_potential)
      library(dplyr)
      library(tidyr)
      table_tree_potential2 <- 
        table_tree_potential %>% 
        dplyr::mutate(pres = ifelse(label %in% colnames(comm_obs), "pres", "abs")) %>% 
        dplyr::mutate(ancestor = table_tree_potential$label[table_tree_potential$parent]) %>% 
        dplyr::mutate(descendant = table_tree_potential$label[table_tree_potential$node]) %>% 
        dplyr::mutate(partition1 = gsub(pattern = ".*_", replacement = "", x = ancestor)) %>% 
        dplyr::mutate(partition2 = gsub(pattern = ".*_", replacement = "", x = descendant)) %>% 
        dplyr::mutate(partition.IS = ifelse(partition1 == "IS" | partition2 == "IS", "IS", NA)) %>% 
        dplyr::mutate(partition.IM = ifelse(partition1 == "IM" & pres == "pres" | partition1 == "ESD" & pres == "pres", "IM", NA)) %>% 
        dplyr::mutate(partition.EM = ifelse(partition1 == "EM" & partition2 == "ESD", "EM", NA)) %>% 
        dplyr::mutate(partition.ESD = ifelse(is.na(partition.IS) == TRUE & is.na(partition.IM) == TRUE & is.na(partition.EM) == TRUE, "ESD", NA))
      
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
      
      PDtotal <- PDinsitu + PDimmigration + PDemmigration + PDexsitu
      
      #joining all results
      data_res <- 
        data.frame(partition = c("PDinsitu", "PDimmigration", "PDemmigration", "PDexsitu", "PDtotal"), 
                   value = c(PDinsitu, PDimmigration, PDemmigration, PDexsitu, PDtotal), community = comm_names[i])
      
      list_res[[i]] <- data_res
    }
    df_res <- do.call(rbind, list_res)
    return(df_res)
  }