# phylogenetic tree
set.seed(42)
phylo <- ape::rcoal(n = 10)
phylo <- ape::makeNodeLabel(phy = phylo)

# community composition matrix
comm <- 
  matrix(sample(c(0, 1), size = 10*20, replace = TRUE),
         nrow = 20, 
         ncol = 10, 
         dimnames = list(paste("comm", 1:20), phylo$tip.label)
  )

# coordinates - this is necessary for "disperal_assembly" algorithm
xy_coords <- 
  matrix(runif(1:10), 
         nrow = nrow(comm),
         ncol = 2, 
         dimnames = list(rownames(comm), c("lng", "lat")
         )
  )

# running sbears with "single_site" algorithm
out_sbears <- 
  calc_sbears(x = comm,
              phy = phylo, 
              coords = xy_coords, 
              method = "single_site")

# sbears object can be used directly in PD_decomposition function

out_decomp <- 
  PD_decomposition(comm = comm,
                   sbears.obj = out_sbears, 
                   phy = phylo, 
                   threshold = 0.5)

# Data frame with results from PD decomposition
out_pd_decomposition <- out_decomp$decomp_potential$PD_decomposition

list_comm_decomposition <- out_decomp$decomp_potential$tree_table_potential
tbl_tree_potential_comm10 <- list_comm_decomposition[[10]]

list_comm_decomposition[[10]] |> View()

library(tidytree)
library(ggtree)

quartz()
ggtree(phy_potential_comm10) + 
  geom_nodelab(aes(label = label)) +
  geom_tiplab()



# Convert tbl_tree back to phylo
list_res_test <- vector(mode = "list", length = length(list_comm_decomposition))
for(i in 1:length(list_comm_decomposition)){
  tbl_tree_potential_comm10 <- list_comm_decomposition[[i]]
  
  tbl_tree_potential_comm10 <- x
  phy_potential_comm10 <- as.phylo(tbl_tree_potential_comm10)
  
  node_times <- branching.times(phy_potential_comm10)
  df_node_times <- 
    data.frame(nodes_names = names(node_times),
               nodes_times = node_times, row.names = NULL)
  
  tbl_tree_potential_comm10_2 <- 
    tbl_tree_potential_comm10 |> 
    left_join(df_node_times, 
              by = c(ancestor = "nodes_names"))
  
  tbl_tree_potential_comm10_3 <- 
    tbl_tree_potential_comm10_2[order(tbl_tree_potential_comm10_2$nodes_times, 
                                      decreasing = TRUE),
    ]
  
  
  tbl_tree_potential_comm10_3_unique <- 
    tbl_tree_potential_comm10_3 |> 
    distinct(ancestor, .keep_all = TRUE)
  
  node_accumulation_potential <- 1:length(tbl_tree_potential_comm10_3_unique$ancestor)
  
  tbl_tree_potential_comm10_4 <- 
    data_frame(tbl_tree_potential_comm10_3_unique, 
               node_accumulation_potential = node_accumulation_potential)
  
  node_accumulation_insitu <- 
    which(tbl_tree_potential_comm10_4$ancestor1 != "IM" & 
            tbl_tree_potential_comm10_4$ancestor1 != "ESD")
  
  tbl_tree_potential_comm10_5 <- 
    tbl_tree_potential_comm10_4 |> 
    mutate(accumulation_insitu = ifelse(ancestor1 != "IM" & ancestor1 != "ESD", TRUE, FALSE)) |> 
    mutate(node_accumulation_model_rec = cumsum(ifelse(tbl_tree_potential_comm10_5$accumulation_insitu, 
                                                       1,
                                                       -1)
    )
    )
  list_res_test[[i]] <- tbl_tree_potential_comm10_5 
}

test_accumm_ltt_comm <- 
  lapply(list_comm_decomposition, function(x){
  # x <- list_comm_decomposition[[1]]
  tbl_tree_potential_comm10 <- x
  phy_potential_comm10 <- as.phylo(tbl_tree_potential_comm10)
  
  node_times <- branching.times(phy_potential_comm10)
  df_node_times <- 
    data.frame(nodes_names = names(node_times),
               nodes_times = node_times, row.names = NULL)
  
  tbl_tree_potential_comm10_2 <- 
    tbl_tree_potential_comm10 |> 
    left_join(df_node_times, 
              by = c(ancestor = "nodes_names"))
  
  tbl_tree_potential_comm10_3 <- 
    tbl_tree_potential_comm10_2[order(tbl_tree_potential_comm10_2$nodes_times, 
                                      decreasing = TRUE),
    ]
  
  
  tbl_tree_potential_comm10_3_unique <- 
    tbl_tree_potential_comm10_3 |> 
    distinct(ancestor, .keep_all = TRUE)
  
  node_accumulation_potential <- 1:length(tbl_tree_potential_comm10_3_unique$ancestor)
  
  tbl_tree_potential_comm10_4 <- 
    data_frame(tbl_tree_potential_comm10_3_unique, 
               node_accumulation_potential = node_accumulation_potential)
  
  node_accumulation_insitu <- 
    which(tbl_tree_potential_comm10_4$ancestor1 != "IM" & 
            tbl_tree_potential_comm10_4$ancestor1 != "ESD")
  
  tbl_tree_potential_comm10_5 <- 
    tbl_tree_potential_comm10_4 |> 
    mutate(accumulation_insitu = ifelse(ancestor1 != "IM" & ancestor1 != "ESD", TRUE, FALSE)) |> 
    mutate(node_accumulation_model_rec = cumsum(ifelse(accumulation_insitu, 
                                                       1,
                                                       -1)
    )
    ) |> 
    mutate(node_accumulation_model_rec2 = ifelse(node_accumulation_model_rec <= 0, 
                                                 0, 
                                                 node_accumulation_model_rec))
  
})

plot(-test_accumm_ltt_comm[[1]]$nodes_times, 
     test_accumm_ltt_comm[[1]]$node_accumulation_model_rec2, type = "line")

for(i in 1:length(test_accumm_ltt_comm)){
  lines(-test_accumm_ltt_comm[[i + 1]]$nodes_times, 
       test_accumm_ltt_comm[[i + 1]]$node_accumulation_model_rec2)
}

lapply(test_accumm_ltt_comm, function(x){
  x$node_accumulation_model_rec2
})


phy_potential_comm10 <- as.phylo(tbl_tree_potential_comm10)

node_times <- branching.times(phy_potential_comm10)
df_node_times <- 
  data.frame(nodes_names = names(node_times),
             nodes_times = node_times, row.names = NULL)

tbl_tree_potential_comm10_2 <- 
  tbl_tree_potential_comm10 |> 
  left_join(df_node_times, 
            by = c(ancestor = "nodes_names"))

tbl_tree_potential_comm10_3 <- 
  tbl_tree_potential_comm10_2[order(tbl_tree_potential_comm10_2$nodes_times, 
                                    decreasing = TRUE),
  ]


tbl_tree_potential_comm10_3_unique <- 
  tbl_tree_potential_comm10_3 |> 
  distinct(ancestor, .keep_all = TRUE)

node_accumulation_potential <- 1:length(tbl_tree_potential_comm10_3_unique$ancestor)

tbl_tree_potential_comm10_4 <- 
  data_frame(tbl_tree_potential_comm10_3_unique, 
           node_accumulation_potential = node_accumulation_potential)

node_accumulation_insitu <- 
  which(tbl_tree_potential_comm10_4$ancestor1 != "IM" & 
          tbl_tree_potential_comm10_4$ancestor1 != "ESD")

tbl_tree_potential_comm10_5 <- 
  tbl_tree_potential_comm10_4 |> 
  mutate(accumulation_insitu = ifelse(ancestor1 != "IM" & ancestor1 != "ESD", TRUE, FALSE)) |> 
  mutate(node_accumulation_model_rec = cumsum(ifelse(tbl_tree_potential_comm10_5$accumulation_insitu, 
                                                     1,
                                                     -1)
                                              )
         )


plot(-tbl_tree_potential_comm10_5$nodes_times, 
     tbl_tree_potential_comm10_5$node_accumulation_model_rec, type = "lines")

comm = comm
sbears.obj = out_sbears 
phy = phylo
threshold = 0.5
make.node.label = TRUE

function(comm, 
         sbears.obj, 
         phy, 
         threshold, 
         make.node.label = TRUE){
  comm <- ifelse(comm >= 1, 1, 0)
  if(make.node.label == TRUE){
    phy <- ape::makeNodeLabel(phy = phy, method = "number", prefix = "Node")
  }
  comm_names <- rownames(comm)
  list_res <- vector(mode = "list", length = nrow(comm)) # object to receive potential tree partition
  list_res2 <- vector(mode = "list", length = nrow(comm)) # object to receive potential tree partition
  list_res3 <- vector(mode = "list", length = nrow(comm)) # object to receive Faith partition
  list_res4 <- vector(mode = "list", length = nrow(comm)) # object to receive Faith partition
  names(list_res) <- rownames(comm)
  
  # reconstruction <- phyloregion::dense2long(t(ifelse(sbears.obj$reconstruction >= threshold, 1, 0))) # nodes predicted from reconstruction
  # reshaping dense to long matrix
  pres_reconstruction <- ifelse(sbears.obj$reconstruction >= threshold, 1, 0)
  reconstruction_df <- as.data.frame(as.table(t(pres_reconstruction)))
  reconstruction <- 
    reconstruction_df[reconstruction_df$Freq >= 1, ]
  reconstruction <- reconstruction[, -3]
  colnames(reconstruction) <- c("grids", "species")
  
  # community <- phyloregion::dense2long(sbears.obj$phylogeny) # Nodes extracted from community phylogeny 
  community <- sbears.obj$site_node_composition # dense matrix
  community <- subset(as.data.frame(as.table(community)), 
                      Freq == 1)[, -3]
  colnames(community) <- c("grids", "species")
  phy_tibble <- tidytree::as_tibble(phy)
  pb <- utils::txtProgressBar(min = 0,      
                              max = length(list_res), 
                              style = 3,    
                              width = 50,   
                              char = "=")   
  for(i in 1:length(list_res)){
    # processing results from reconstruction
    # i = 20
    nodes_reconstruction <-
      reconstruction |> 
      subset(grids == comm_names[i]) # model reconstruction
    
    nodes_comm <- 
      community |> 
      subset(grids == comm_names[i]) # from phylogeny based on current species occupancy
    spp_comm <- names(which(comm[i, ] == 1)) # species in community
    comm_tible <- 
      phy_tibble[match(c(spp_comm, as.character(nodes_comm$species)), 
                       phy_tibble$label),
      ]
    offspring_data <-
      tidytree::offspring(phy_tibble, 
                          as.character(nodes_reconstruction$species))
    if(class(offspring_data)[1] == "list"){
      df_offspring_rec <- 
        do.call(rbind, 
                tidytree::offspring(phy_tibble, 
                                    as.character(nodes_reconstruction$species)
                )
        )
    } else{
      df_offspring_rec <- offspring_data
    }
    
    # starting calculation of PD components
    if(is.null(df_offspring_rec) == TRUE){ # when there is no ancestor occupuying communities
      tree_potential <- ape::keep.tip(phy = phy, tip = spp_comm)
      spp_potential <- spp_comm
      spp_potential_all <- spp_comm
    } else{
      spp_potential <- 
        df_offspring_rec |> 
        dplyr::distinct(label) # species in observed communities
      spp_potential2 <- phy$tip.label[phy$tip.label %in% spp_potential$label] # species names from nodes estimated in reconstruction
      spp_potential_all <- unique(c(spp_potential2, spp_comm)) # joining community and reconstruction 
      tree_potential <- ape::keep.tip(phy = phy, tip = spp_potential_all) # keeping only species observed in communities and estimated in reconstruction
    }
    
    # naming node categories
    # ISDiv - In Situ Diversification nodes
    nodes_insitu <- intersect(as.character(nodes_reconstruction$species), 
                              as.character(nodes_comm$species))
    tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)] <-
      paste(tree_potential$node.label[match(nodes_insitu, 
                                            tree_potential$node.label)], 
            "IS", sep = "_") # naming as IS standing for In Situ
    
    # IM - Immigration nodes
    node_immigration <- setdiff(as.character(nodes_comm$species),
                                as.character(nodes_reconstruction$species))
    
    # Tagging nodes with IM - Immigration nodes
    tree_potential$node.label[match(node_immigration, tree_potential$node.label)] <-
      paste(tree_potential$node.label[match(node_immigration, 
                                            tree_potential$node.label)], "IM", 
            sep = "_") 
    
    # EM - emigration nodes
    node_emigration <- setdiff(as.character(nodes_reconstruction$species),
                               as.character(nodes_comm$species))
    tree_potential$node.label[match(node_emigration, tree_potential$node.label)] <-
      paste(tree_potential$node.label[match(node_emigration, 
                                            tree_potential$node.label)], "EM", sep = "_")
    
    # ESDiv - Ex Situ Diversification nodes
    node_exsitu <- tree_potential$node.label[-grep(pattern = "_", tree_potential$node.label)] 
    tree_potential$node.label[match(node_exsitu, tree_potential$node.label)] <-
      paste(tree_potential$node.label[match(node_exsitu, 
                                            tree_potential$node.label)], 
            "ESD", sep = "_")
    
    # Phylogenetic table tree with potential lineages in communities
    table_tree_potential <- tidytree::as_tibble(tree_potential)
    
    # Phylogenetic tree with observed lineages in communities
    tree_faith <- ape::keep.tip(phy = tree_potential, tip = spp_comm) # keeping only species observed in communities and estimated in reconstruction
    
    # Phylogenetic table tree for observed lineages in communities
    table_tree_faith <- tidytree::as_tibble(tree_faith)
    
    # Community matrix that will be used to calculate PD for observed lineages
    comm_faith <- 
      matrix(rep(1, length(spp_comm)), nrow = 1, ncol = length(spp_comm), 
             dimnames = list(comm_names[i], spp_comm))
    
    # Community matrix that will be used to calculate potential PD in communities
    comm_potential <- 
      matrix(rep(1, length(spp_potential_all)), nrow = 1, ncol = length(spp_potential_all),
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
    utils::setTxtProgressBar(pb, i)
  }
  
  # summarizing results
  decomp_potential <- organize_pd_decomp(list1.decomp = list_res, list2.treetable = list_res2) # function used to organize the results from pd decomposition
  decomp_faith <- organize_pd_decomp(list1.decomp = list_res3, list2.treetable = list_res4)
  list_res_final <- vector(mode = "list")
  list_res_final$decomp_potential <- decomp_potential
  list_res_final$decomp_faith <- decomp_faith
  return(list_res_final)
}