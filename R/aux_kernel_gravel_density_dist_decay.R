#' Occupancy probability of each node in each site considering a kernel density funcion
#'
#' @param x Node by communities matrix
#' @param phy phylo object 
#' @param coords data frame with spatial coordinate of communities
#'
#' @returns A list with length equal the number of sites. Each element of the list
#'     contains a matrix with nodes in rows and sites in columns and the probability
#'     of occurrence of each ancestral node in each site.  
#' @export

comp_kernel_gravel <- 
  function(x, phy, coords, w_slope, min_disp_prob){
    
    # modified site x node matrix to receive dispersal assembly results
    node.anc.area.spat <- matrix(NA, 
                                 nrow = phy$Nnode, 
                                 ncol = ncol(x),
                                 dimnames = list(phy$node.label, colnames(x)))
    
    # pairwise geographical distance matrix in km
    r <- 
      scales::rescale(geodist::geodist(x = coords, measure = "geodesic")/
                        1000, diag = T, upper = T, c(0, 1)) # in km
    rownames(r) <- colnames(r) <- colnames(x) # naming distance matrix with assemblage names
    max_disp_dist <- sqrt(-log(min_disp_prob)/w_slope) # max distance to be used in kernel 
    anc_list <- vector(mode = "list", length = nrow(r)) # list with all distance decay values for each focal site
    site_values_pernode <- matrix(NA, 
                                  ncol(x), 
                                  ncol(x),
                                  dimnames = list(colnames(x), colnames(x)))
     # matrix to receive distance decay values of each node per focal cell
    site_values_pernode <- 
      matrix(NA, 
             nrow = ncol(x),
             ncol = ncol(x),
             dimnames = list(colnames(x), colnames(x)))
    
    # distance decay for all focal cells for all nodes each element of the list is a focal cell
    for (i in 1:nrow(r)){
      r_below_threshold <- which(r[i,] <= max_disp_dist)
      r_pruned <- r[r_below_threshold, r_below_threshold]
      dist.decay <- matrix(NA, 
                           nrow = phy$Nnode, 
                           ncol = nrow(r_pruned),
                           dimnames = list(phy$node.label, 
                                           rownames(r_pruned)
                           )
      )
      
      for (k in rownames(r_pruned)){
        for (j in 1:phy$Nnode){ # each node in all sites
          for (p in colnames(r_pruned)){
            dist.decay[j, p] <- x[j, k]*exp(1)^-(w_slope*r_pruned[k, p]^2)
          }
        }
      }
      anc_list[[i]] <- dist.decay
    }
    
    # 
    for (j in 1:phy$Nnode){ # node probability of occupancy in each cell from reconstruction
      # j = 1
      for(l in 1:length(anc_list)){ # distance decay for all nodes from a focal cell to all other cells
        dist.decay <- anc_list[[l]] 
        for (k in colnames(dist.decay)){
          site_values_pernode[l, k] <- anc_list[[l]][j, k]
        }
      }
      for(k in 1:ncol(x)){
        site_values_pernode_site <- site_values_pernode[,k][which(!is.na(site_values_pernode[,k]==TRUE))]
        dens_site_node_pernode <- density(site_values_pernode_site, from = 0, to = 1)
        #dens_site_node_pernode <- 
        #  density(site_values_pernode[, k],
        #          from = 0, 
        #          to = 1, 
        #          na.rm = TRUE)
        node.anc.area.spat[j, k] <- 
          dens_site_node_pernode$x[which.max(dens_site_node_pernode$y)]
      }
    }
    res <- vector(mode = "list", length = 2)
    res$list_focal_cell_dist_decay <- anc_list # distance decay object
    res$node_anc_area_spat <- node.anc.area.spat # probability of dispersion 
    return(res)
  }    


