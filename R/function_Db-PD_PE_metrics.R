#' Diversification-based PD and PE
#' 
#' @details This function computes two modified versions of Phylogenetic Diversity (PD) and Phylogenetic
#'     endemism by considering in the calculation the in-situ diversification. the Db-PD and Db-PE are calculated
#'     by coupling in the calculation an ancetral area reconstruction that divide the phylogenetic tree in two parts, 
#'     one that correspond to dispersion events and another that correspond to in-situ diversification events. 
#'     We use the in-situ diversification part to calculate both Db-PD and Db-PE
#' 
#' @param W Occurrence matrix, rows are assemblages and columns are species
#' @param tree Phylogenetic hipothesis in newick format
#' @param ancestral.area One column data frame with nodes in rows and one column indicating the occurrence area of nodes
#' @param biogeo One columns data frame with rows indicating assemblages and one column indicating the biome/ecoregion of each assemblage
#' @param PD Logical, if TRUE (default) Db-PD will be computed
#' @param PE Logical, if TRUE (default) Db-PE will be computed
#'
#' @return A data frame containing the values of original PD and PE and also their model
#'     based version
#' 
#' @author Gabriel Nakamura <gabriel.nakamura.souza@@gmail.com>
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' W_toy<- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0),
#' nrow= 3,
#' ncol= 5,
#' dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
#' c(paste("s", 1:5, sep=""))))
#'  
#'  biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
#'  ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
#'  assemblage_phylo_metrics <- calc_div_based_metrics(W_toy, toy_treeEx, ancestral_area_toy, biogeo_toy)
#' 
#' }
calc_div_based_metrics <- function(W,
                              tree,
                              ancestral.area, 
                              biogeo,
                              PD = TRUE,
                              PE = TRUE){
  
  # total phylogenetic diversity and endemism -------------------------------
  if(!is.matrix(W) == TRUE){
    if(is.data.frame(W) == TRUE){
      W <- as.matrix(W)
      rownames(W) <- 1:nrow(W)
    } else{
      stop("W must be a occurrence matrix with presences (1) or absences (0)")
    }
  }
  PDt <- picante::pd(samp = W, tree = tree)$PD # faster option
  PEt <- phyloregion::phylo_endemism(x = W , phy = tree, weighted = TRUE) #Phylogenetic endemism sensu Rosauer
  
  # basic information from nodes and ancestral reconstruction ---------------
  nodes.list <- get_nodes_info_core(W = W, tree = tree, ancestral.area = ancestral.area, biogeo = biogeo)
  
  names_spp_noNull<- lapply(
    lapply(
      lapply(nodes.list,
             function(x){
               is.na(x$nodes_species)
             }
      ),
      function(x){
        which(x == FALSE)
      }
    ),
    names)
  nodes_species_noNull <- vector(mode = "list", length = nrow(W))
  for(i in 1:length(names_spp_noNull)){
    nodes_species_noNull[[i]] <- nodes.list[[i]]$nodes_species[names_spp_noNull[[i]]]
  }
  
  ####organize node matrix
  nodes_species_noNull_org <- nodes_species_noNull
  list_matrix_nodes <- vector(mode = "list", length = length(nodes_species_noNull_org))
  for(i in 1:length(nodes_species_noNull)){
    #i= 488
    if(length(nodes_species_noNull[[i]]) == 0){
      list_matrix_nodes[[i]] <- NA
    } else{
      names_spp <- names(nodes_species_noNull[[i]])
      list_nodes_org <- vector(mode = "list", length = length(names_spp))
      for(j in 1:length(names_spp)){
        #j= 1
        matrix_nodesSpp_nonull <- matrix(NA, nrow = round(length(unlist(nodes_species_noNull[[i]][j]))), ncol= 2)
        nodes_org <- c(sort(nodes_species_noNull[[i]][j][[1]],
                            decreasing = FALSE),
                       which(tree$tip.label == names_spp[j]))
        for(k in 1:nrow(matrix_nodesSpp_nonull)){
          matrix_nodesSpp_nonull[k, ] <- c(nodes_org[k], nodes_org[k + 1])
        }
        list_nodes_org[[j]]<- matrix_nodesSpp_nonull
      }
      list_matrix_nodes[[i]]<- list_nodes_org
    }
  }
  
  
  
  # extracting all internal insitu branch length ----------------------------
  
  list_matrix_edges <- vector(mode = "list", length= length(list_matrix_nodes))
  for(i in 1:length(list_matrix_nodes)){
    if(any(is.na(list_matrix_nodes[[i]])) == TRUE){
      list_matrix_edges[[i]] <- NA
    } else{
      list_matrix_edges[[i]] <- unique(do.call(rbind, 
                                               list_matrix_nodes[[i]]))
    }
  }
  
  
  list_brLength_divLocal <- vector(mode = "list", length = length(list_matrix_edges))
  for(i in 1:length(list_matrix_edges)){
    if(any(is.na(list_matrix_edges[[i]])) == TRUE){
      list_brLength_divLocal[[i]] <- NA
    } else{
      pos_insitu_nodes<- apply(tree$edge, MARGIN = 1, function(x){
        which(x == apply(as.matrix(list_matrix_edges[[i]]), MARGIN = 1, function(l) l))
      })
      list_brLength_divLocal[[i]]<- tree$edge.length[which(unlist(lapply(pos_insitu_nodes,
                                                                         function(x) length(x) > 1)
      ) == TRUE)
      ]
    }
  }
  
  
  # computing diversification-based PD --------------------------------------
  
  if(PD == TRUE){
    
    PDlocal <- unlist(lapply(list_brLength_divLocal,
                             function(x){
                               sum(x)
                             }
    )
    )
    names(PDlocal)<- rownames(W) #PDlocal
  }
  
  
  
  # computing diversification-based PE --------------------------------------
  
  if(PE == TRUE){
    names(list_matrix_edges) <- paste("site", 1:nrow(W), sep= "_")
    all_comb <- gtools::permutations(length(list_matrix_edges), 2,
                                     names(list_matrix_edges))
    list_occurence<- vector(mode= "list", length= length(list_matrix_edges))
    
    #calculating the number of branches co-occurence
    list_matrix_edges_allTips <- vector(mode = "list", length = length(list_matrix_edges))
    names(list_matrix_edges_allTips)<- names(list_matrix_edges)
    for(i in 1:length(list_matrix_edges)){
      #i= 42
      tip_edges<- do.call(rbind, lapply(ape::nodepath(tree)[which(W[i,] == 1)], function(x){
        x[(length(x) - 1):length(x)]
      }))
      if(is.na(list_matrix_edges[[i]])){
        list_matrix_edges_allTips[[i]]<- unique(tip_edges) # only tip edges
      } else{
        list_matrix_edges_allTips[[i]]<- unique(rbind(list_matrix_edges[[i]], tip_edges)) # tip and internal branches
      }
    }
    
    for(i in 1:length(list_matrix_edges_allTips)){
      posit<- which(all_comb[, 1] == names(list_matrix_edges_allTips)[i])
      comb_occ<- matrix(NA, nrow = nrow(list_matrix_edges_allTips[[i]]), ncol = length(list_matrix_edges_allTips)-1)
      for(j in 1:length(posit)){
        if (any(is.na(list_matrix_edges_allTips[[all_comb[posit[j], 2]]])) == TRUE){
          comb_occ[, j]<- FALSE
        } else {
          comb_occ[,j]<- list_matrix_edges_allTips[[all_comb[posit[j], 1]]][, 1] %in% list_matrix_edges_allTips[[all_comb[posit[j], 2]]][, 1] &
            list_matrix_edges_allTips[[all_comb[posit[j], 1]]][, 2] %in% list_matrix_edges_allTips[[all_comb[posit[j], 2]]][, 2]
        }
      }
      list_occurence[[i]]<- comb_occ
    }
    
    
    # naming list_occurrence
    for(i in 1:length(list_occurence)){
      if(any(is.na(list_occurence[[i]]))){
        list_occurence[[i]]<- NA
      } else{
        colnames(list_occurence[[i]]) <- biogeo[-i, 1]
      }
    }
    
    # binding occurrence and branch length
    list_edge_brLen_insitu <- vector(mode = "list", length = length(list_occurence))
    for(i in 1:length(list_occurence)){
      list_edge_brLen_insitu[[i]]<- cbind(list_brLength_divLocal[[i]], list_occurence[[i]])
    }
    
    ####correcting denominator of PE for ancestral state occurrence
    ancestral_nodes <- cbind(ancestral.area, node = rownames(ancestral.area)) # binding nodes name with ancestral state
    
    
    ### ancestral state for each node that compose each community
    ancestral_comm_temp<- lapply(
      lapply(list_matrix_edges_allTips,
             function(x){
               if(any(is.null(dim(x))) == TRUE){
                 NA
               } else {
                 apply(x, MARGIN = 1, function(l){
                   cbind(ancestral_nodes[match(l[1], ancestral_nodes[, 2]), 1], ancestral_nodes[match(l[2], ancestral_nodes[, 2]), 1])
                 }
                 )
               }
             }
      ),
      function(m) t(m)
    )
    
    ### binding nodes with their respective ancestral state
    list_matrix_edges_AS<- lapply(list_matrix_edges_allTips, function(x){
      if(any(is.na(x)) == TRUE){
        NA
      } else{
        cbind(x, rep(NA, nrow(x)), rep(NA, nrow(x)))
      }
    }
    )
    for(i in 1:length(list_matrix_edges_AS)){
      if(any(is.null(dim(list_matrix_edges_AS[[i]]))) == TRUE){
        list_matrix_edges_AS[[i]]<- NA
      } else{
        list_matrix_edges_AS[[i]][, c(3,4)]<- ancestral_comm_temp[[i]]
      }
    }
    
    # adding ancestral occurrence and branch length
    list_all_DbInfo<- vector(mode = "list", length= length(list_brLength_divLocal))
    for(i in 1:length(list_matrix_edges_AS)){
      if(any(is.na(list_brLength_divLocal[[i]])) == TRUE){
        list_edge_brLen_insitu[[i]][, 1] <- tree$edge.length[tree$edge[, 1] %in%  list_matrix_edges_allTips[[i]][, 1] &
                                                               tree$edge[, 2] %in%  list_matrix_edges_allTips[[i]][, 2]
        ] # obtaining branch lengths for terminal tips
        denom_tmp2<- rowSums(list_occurence[[i]]) # denominator for terminal tips
        list_all_DbInfo[[i]]<- cbind(list_edge_brLen_insitu[[i]][,1], ((denom_tmp2) + 1)) # joining information for terminal tips
      } else{ # correcting the denominator for occurrence of internal branch lengths
        a<- strsplit(list_matrix_edges_AS[[i]][, 3], split= "")
        b<- strsplit(list_matrix_edges_AS[[i]][, 4], split= "")
        denom_tmp2<- numeric()
        for(j in 1:length(a)){
          if(any(is.na(b[[j]])) == TRUE){
            denom_tmp2[j] <- sum(list_occurence[[i]][j,]) # denominator for tip edges
          } else{
            denom_tmp2[j]<- sum(list_edge_brLen_insitu[[i]][j, ][!is.na(match(names(list_edge_brLen_insitu[[i]][j, ]),
                                                                              a[[j]][which(a[[j]] == b[[j]])]
            )
            )
            ]
            ) # denominator that counts only for bioregions that match with current occurrence
          }
        }
      }
      list_all_DbInfo[[i]] <- cbind(list_edge_brLen_insitu[[i]][,1], ((denom_tmp2) + 1)) # binding denominator and branch length for each edge
    }
    
    #### Calculating Db_PE
    Db_PE_res <- lapply(
      lapply(list_all_DbInfo,
             function(x){
               if(any(is.null(dim(x))) == TRUE){
                 NA
               } else{
                 apply(x, MARGIN = 1,
                       function(l){
                         if(any(is.na(l)) == T){
                           NA
                         } else{
                           l[1]/l[2]
                         }
                       }
                 )
               }
             }
      ),
      sum)
    res_Db_PE <- matrix(unlist(Db_PE_res), nrow= length(list_matrix_edges_AS), ncol= 1,
                        dimnames = list(names(list_matrix_edges_AS), "Db_PE"))
  }
  
  # Output for community metrics ---------------------------------------
  
  if(PD == TRUE & PE == TRUE){
    Db_comm_metrics <- data.frame(PE = PEt, Db_PE = res_Db_PE, PD = PDt, PDlocal= PDlocal)
  }
  if(PD == TRUE & PE == FALSE){
    Db_comm_metrics <- data.frame(PD = PDt, PDlocal = PDlocal)
  }
  if(PE == TRUE & PD == FALSE){
    Db_comm_metrics <- data.frame(PE= PEt, Db_PE = res_Db_PE)
  }
  
  return(Db_comm_metrics)
  
}

