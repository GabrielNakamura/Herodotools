# primates
res_primates <- load(here::here("inst", "extdata", "primates_PD_decomp.RData")) # results from primates
res_ada_primates <- load(here::here("inst", "extdata", "Primates_ada_test2.RData"))
comm = primates.PAM2 # community matrix
ada.obj = ada.obj2 # output from ada reconstruction
phy = primate.tree # phylogenetic tree
threshold = 0.9  # threshold to be used in function

# tiranideos

res_tiranideo <- load(here::here("inst", "extdata", "ada_tiranideos.RData")) # results from primates
comm = comm # community matrix
ada.obj = ada.obj # output from ada reconstruction
phy = phy # phylogenetic tree
threshold = 0.9  # threshold to be used in function


PD_decomposition <- 
  function(comm, ada.obj, phy, threshold){
    comm <- ifelse(comm >= 1, 1, 0)
    phy <- ape::makeNodeLabel(phy = phy, method = "number", prefix = "Node")
    comm_names <- rownames(comm)
    list_res <- vector(mode = "list", length = nrow(comm))
    list_res2 <- vector(mode = "list", length = nrow(comm))
    names(list_res) <- rownames(comm)
    
    reconstruction <- phyloregion::dense2long(t(ifelse(ada.obj$reconstruction >= threshold, 1, 0))) # nodes predicted from reconstruction
    community <- phyloregion::dense2long(ada.obj$phylogeny) # Nodes extracted from community phylogeny 
    phy_tibble <- tidytree::as_tibble(phy)
    for(i in 2001:3100){
      # i = 500
      # setting a progress bar
      pb <- txtProgressBar(min = 0,      
                           max = length(list_res), 
                           style = 3,    
                           width = 50,   
                           char = "=")   
      
      # processing results from reconstruction
      nodes_reconstruction <-
        reconstruction %>% 
        subset(grids == comm_names[i]) # reconstruction
      nodes_comm <- 
        community %>% 
        subset(grids == comm_names[i]) # phylogeny
      
      # finding the most basal node among all set of nodes including reconstruction and community nodes for a given community
      node_sequence <- gsub(pattern = "Node", replacement = "", unique(c(nodes_reconstruction$species, nodes_comm$species)))
      data_nodes <- data.frame(nodes = unique(c(nodes_reconstruction$species, nodes_comm$species)), node_sequence = as.numeric(node_sequence))
      nodes_all <- phy_tibble[phy_tibble$label %in% data_nodes$nodes, "node"]$node
      spp_potential <- unique(unlist(lapply(nodes_all, function(x) phytools::getDescendants(tree = phy, node = x)))) %in% 1:length(phy$tip.label)
      desc_potential <- unique(unlist(lapply(nodes_all, function(x) phytools::getDescendants(tree = phy, node = x))))[spp_potential]
      mrca_potential <- ape::getMRCA(phy = phy, tip = desc_potential)
      tree_potential <- ape::extract.clade(phy = phy, node = mrca_potential)
      
      
      # naming node categories
      # ISDiv - In Situ Diversification nodes
      nodes_insitu <- intersect(nodes_reconstruction$species, nodes_comm$species)
      tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)], "IS", sep = "_") # nomeando nos com IS
      
      # IM - Immigration nodes
      node_immigration <- setdiff(nodes_comm$species, nodes_reconstruction$species)
      tree_potential$node.label[match(node_immigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_immigration, tree_potential$node.label)], "IM", sep = "_")
      
      # EM - Emmigration nodes
      node_emmigration <- setdiff(nodes_reconstruction$species, nodes_comm$species)
      tree_potential$node.label[match(node_emmigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_emmigration, tree_potential$node.label)], "EM", sep = "_")
      
      # ESDiv - Ex Situ Diversification nodes
      node_exsitu <- tree_potential$node.label[-grep(pattern = "_", tree_potential$node.label)] 
      tree_potential$node.label[match(node_exsitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_exsitu, tree_potential$node.label)], "ESD", sep = "_")
      
      table_tree_potential <- tidytree::as_tibble(tree_potential)
      comm_obs <- names(which(comm[comm_names[i], ] == 1))
      # Organizing data to calculate PD components
      table_tree_potential2 <- 
        table_tree_potential %>% 
        dplyr::mutate(pres = ifelse(label %in% comm_obs, "pres", "abs")) %>% 
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
                                              ancestor1 == "IM" & descendent1 == "EM" |
                                              ancestor1 == "ESD" & descendent1 == "IS" |
                                              ancestor1 == "ESD" & descendent1 == "EM",
                                            "IM", NA)) %>% 
        dplyr::mutate(partition.EM = ifelse(ancestor1 == "IS" & descendent1 == "IM" |
                                              ancestor1 == "EM" & descendent1 == "ESD" |
                                              ancestor1 == "EM" & descendent1 == "IM" | 
                                              ancestor1 == "IS" & descendent1 == "ESD" |
                                              ancestor1 == "EM" & pres == "abs" & descendent1 != "IS", 
                                            "EM", NA)) %>% 
        dplyr::mutate(partition.ESD = ifelse(ancestor1 == "ESD" & descendent1 == "ESD" |
                                               ancestor1 == "ESD" & descendent1 == "IM" |
                                               ancestor1 == "ESD" & pres == "abs", "ESD", NA)) %>% 
        dplyr::mutate(partition.unknown = ifelse(is.na(partition.IS) & is.na(partition.IM) & is.na(partition.EM) & is.na(partition.ESD), "unknown", "known")) %>% 
        dplyr::mutate(class.unknown = dplyr::case_when(partition.unknown == "unknown" ~ ifelse(pres == "pres", paste(ancestor1, pres, sep = "_"), paste(ancestor1, descendent1, sep = "_")))) %>% 
        dplyr::mutate(class.known = dplyr::case_when(partition.unknown == "known" ~ ifelse(pres == "pres", paste(ancestor1, pres, sep = "_"), paste(ancestor1, descendent1, sep = "_"))))
     
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
      
      # progress bar
      setTxtProgressBar(pb, i)
    }
    
    # summarizing results
    df_res <- do.call(rbind, list_res)
    list_res_final <- vector(mode = "list")
    list_res_final$PD_decomposition <- df_res
    list_res_final$tree_table_potential <- list_res2
    
    # calculating unknown node links for each community
    
    list_unknown_count <- 
      lapply(list_res_final$tree_table_potential[1:10], function(x){
      x %>% 
        filter(partition.unknown == "unknown") %>% 
        group_by(class.unknown) %>% 
        add_count() %>% 
        distinct() %>% 
        select(parent, node, branch.length, label, ancestor, descendant, class.unknown, n)
    })
    list_res_final$unknown <- list_unknown_count
    
    # matrix of PD decomposition - dense format
    
    PD_decomposition_wide <- tidyr::pivot_wider(data = list_res_final$PD_decomposition, names_from = "partition", values_from = "value")
    list_res_final$dense_matrix_PD <- PD_decomposition_wide
    
    return(list_res_final)
  }