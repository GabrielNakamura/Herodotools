
# this is a half length case lineage, the change is at the half of the length

W[c(593, 594, 610, 626, 642), ] # species A juninensis in D site - this species is BD and the next node is CD

W[which(biogeo[,1] == "D"), "A_juninensis"]


l.Jetz.local <- lapply(nodes.list , function(site){
  lapply(1:length(nodes.list[[593]]$nodes_species), function(j){
    # nodes_div <- nodes.list[[3]]$nodes_species[[1]]
    # sp <- names(nodes.list[[3]]$nodes_species)[1]
    # j = 1
    nodes_div <- nodes.list[[593]]$nodes_species[[j]]
    sp <- names(nodes.list[[593]]$nodes_species)[j]
    
    if(length(nodes_div) == 1){
      # Firts testing if the node is from an anagenetic origin
      node_test_anagenetic <- tree$node.label[stringr::str_detect(tree$node.label, as.character(nodes_div))]
      if(grepl("\\bnew_", node_test_anagenetic) == TRUE){
        # if TRUE we get the height of the node up to the first cladogenetic node
        internal.brlen_div <- tree$edge.length[which(tree$edge[,
                                         1] %in% nodes_div)]
      } else{
        #  get the internal branch length from the parent node to the target node
        internal.brlen_div <- tree$edge.length[which(tree$edge[,
                                                               2] %in% nodes_div)] #branch lengths (in times) for internal branch lengths of the most ancient ancestral
      }
     
    } else{
      nodes_div <- nodes_div[-1] # this is used to remove the path to the node in another area - ancient ancestral that was presented at Ecoregion of local i to species j
      internal.brlen_div <- tree$edge.length[which(tree$edge[,
                                                             2] %in% nodes_div)] #ages for internal nodes of nodes_div object
    }
    if (length(internal.brlen_div) != 0) { #starts the calculation for ED measure from Redding et al.
      internal.brlen_div <- internal.brlen_div * switch(type, equal.splits = sort(rep(0.5,
                                                                                      length(internal.brlen_div))^c(1:length(internal.brlen_div))),
                                                        fair.proportion = {
                                                          for (j in 1:length(nodes_div)) {
                                                            sons <- picante::.node.desc(tree, nodes_div[j])
                                                            n.descendents <- length(sons$tips)
                                                            if (j == 1) portion <- n.descendents else portion <- c(n.descendents,
                                                                                                                   portion)
                                                          }
                                                          1/portion
                                                        })
      
    }
    
    if(any(is.na(nodes_div))){
      JetzLocal <- 0.00001
    }else{
      ED_div <- sum(internal.brlen_div,
                    tree$edge.length[ape::which.edge(tree, sp)]) #Local ED - modified ED considering only the edges of ancestors inside the biogeo of local i for species j
      
      EDtotal_spp <- EDtotal$w[which(EDtotal$Species == sp)]
      Jetz_total_spp <- Jetz_total[sp] # Diversification calculated according to Jetz for species j
      JetzLocal <- (Jetz_total_spp*(ED_div/EDtotal_spp)) # Modified Jetz to calculate only for local diversification
      
    }
    JetzLocal #Jetz local diversification
    
  })#lapply j == species
})