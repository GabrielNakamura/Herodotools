library(dplyr)
devtools::load_all()

# primates
res_primates <- load(here::here("inst", "extdata", "primates_PD_decomp.RData")) # results from primates
res_ada_primates <- load(here::here("inst", "extdata", "Primates_ada_test2.RData"))
comm = primates.PAM2 # community matrix
ada.obj = ada.obj2 # output from ada reconstruction
phy = primate.tree # phylogenetic tree
threshold = 0.9  # threshold to be used in function

# tiranideos

res_tiranideo <- load(here::here("inst", "extdata",  "Tyrannidae.RData")) # results from primates
comm = comm # community matrix
ada.obj = ada.obj.tyr # output from ada reconstruction
phy = phy # phylogenetic tree
threshold = 0.9  # threshold to be used in function

# sigmodontideos

res_sigmo <- load(here::here("inst", "extdata", "ada_sigmo", "adaobj_sig.RData")) # results from primates
comm = comm # community matrix
ada.obj = ada.obj.sig # output from ada reconstruction
phy = phy # phylogenetic tree
threshold = 0.9  # threshold to be used in function


# Muroidea
res_muro <- load(here::here("inst", "extdata", "Muroidea_test.RData")) # results from primates
comm = comm # community matrix
comm <- comm[ -which(rowSums(comm) == 0 | rowSums(comm) == 1 | rowSums(comm) == 2), ]
ada.obj = ada.obj.mur # output from ada reconstruction
phy = phy # phylogenetic tree
threshold = 0.9  # threshold to be used in function

test_muro <- PD_decomposition(comm = comm.test, ada.obj = ada.obj.mur, phy = phy, threshold = 0.9)
test2 <- readRDS(here::here("inst", "extdata", "res_muroidea.rds"))
test2$decomp_potential$dense_matrix_PD
test2$decomp_faith$dense_matrix_PD
### Versão leandro PD decomposition

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
    for(i in 1:length(list_res)){
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
                                              ancestor1 == "ESD" & descendent1 == "EM"|
                                              ancestor1 == "ESD" & pres == "pres"| 
                                              ancestor1 == "IM" & pres == "pres",
                                            "IM", NA)) %>% 
        dplyr::mutate(partition.EM = ifelse(ancestor1 == "IS" & descendent1 == "IM" |
                                              ancestor1 == "EM" & descendent1 == "ESD" |
                                              ancestor1 == "EM" & descendent1 == "IM" | 
                                              ancestor1 == "IS" & descendent1 == "ESD" |
                                              ancestor1 == "IS" & pres == "abs" | 
                                              ancestor1 == "EM" & pres == "abs", 
                                            "EM", NA)) %>% 
        dplyr::mutate(partition.ESD = ifelse(ancestor1 == "ESD" & descendent1 == "ESD" |
                                               ancestor1 == "ESD" & descendent1 == "IM" |
                                               ancestor1 == "IM" & descendent1 == "IM" |
                                               ancestor1 == "IM" & descendent1 == "ESD" |
                                               ancestor1 == "IM" & pres == "abs" | 
                                               ancestor1 == "ESD" & pres == "abs", "ESD", NA))
      
      # adding a group to all species - this will be useful to plot the partitions
      table_tree_potential2 <- 
        table_tree_potential2 %>% 
        mutate(group = coalesce(partition.IS, partition.IM, partition.EM, partition.ESD))  
      
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
      
      PDemigration <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.EM == "EM") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = TRUE)
      
      PDexsitu <- 
        table_tree_potential2 %>% 
        dplyr::filter(partition.ESD == "ESD") %>% 
        dplyr::select(branch.length) %>% 
        sum(na.rm = TRUE)
      
      PDtotal <- PDinsitu + PDimmigration + PDemigration + PDexsitu # adding up all components
      PDfaith <- picante::pd(samp = t(as.matrix(comm[i, ])), tree = phy,include.root = FALSE)[1, 1]  # including PD calculated accordingly to Faith approach
      
      #joining all results
      data_res <- 
        data.frame(partition = c("PDinsitu", "PDimmigration", "PDemigration", "PDexsitu", "PDtotal", "PDfaith"), 
                   value = c(PDinsitu, PDimmigration, PDemigration, PDexsitu, PDtotal, PDfaith), community = comm_names[i])
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
    
    # matrix of PD decomposition - dense format
    PD_decomposition_wide <- tidyr::pivot_wider(data = list_res_final$PD_decomposition, names_from = "partition", values_from = "value")
    list_res_final$dense_matrix_PD <- PD_decomposition_wide
    
    return(list_res_final)
  }

#######

comm <- comm[-which(rowSums(comm) == 0), ]

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
    # list_res <- list_res[600:700]
    pb <- txtProgressBar(min = 0,      
                         max = length(list_res), 
                         style = 3,    
                         width = 50,   
                         char = "=")   
    for(i in 1:length(list_res)){
      # i = 600
      # setting a progress bar
      
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