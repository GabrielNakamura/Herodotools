#' Site-Based Estimation of Ancestral Range of Species
#'
#' @param x Matrix object with community/assemblage composition.
#'     Rows represent sites and columns represents species.
#' @param phy A phylogenetic tree of \code{phylo} class object.
#' @param coords A rectangular object (it can be a matrix, data.frame or tibble)
#'     containing two columns, one with longitude and another with latitude 
#'     coordinates in this order.
#' @param method character indicating how ancestral range probabilities are computed.
#'     The options are "single_site" or "disp_assembly"
#' @param w_slope A scalar representing the slope of the dispersal kernel function.
#' @param compute.node.by.sites Logical, TRUE (default) computes a matrix of node occurrence by site.
#' @param make.node.label Logical, if TRUE (default) the nodes of the phylogeny will be named as the letter "N" preceding node number
#' @details
#' SBEARS (Site-Based Estimation of Ancestral Range of Species) is a method for 
#' ancestral state reconstruction of species geographic ranges. It operates at 
#' fine spatial resolution and does not require predefined discrete biogeographic 
#' areas. The method estimates ancestral ranges under a maximum likelihood 
#' framework, using site-level occurrence information.
#' \itemize{
#'     \item the method used to compute ancestral range probabilities can be 
#'         \code{single_site} or \code{disp_assembly}. For both methods SBEARS use the function 
#'         \code{\link[phytools]{fastAnc}} to reconstruct the the ancestral ranges.
#'         In \code{single_site} each site is taken 
#'         as a binary categorical trait representing the presence (state 1) or 
#'         the absence (state 0) of every species in the site, which is reconstructed
#'         along the nodes of the phylogeny.
#'         In \code{disp_assembly} option the single site reconstruction is 
#'         carried out, producing a matrix describing standardized probabilities
#'         of occurrence of each node at each site. The difference here is that 
#'         an additional step is performed to compute the probability of a cell 
#'         being colonized by an ancestral species originating from a focal cell.
#'         The probability of colonization decreases increasing the distance from
#'         the focal cell. This probability of an ancestral species to disperse
#'         from a cell to any other can be defined using a dispersal kernel function
#'         proposed by Gravel et al (2006).
#' }
#' 
#' @return a list with three elements. reconstruction is the result of ancestral
#'.   area reconstruction; phylogeny is the matrix containing
#'    the occurrence of nodes in sites and joint.phylo.obs is the joint occurrence 
#'    of nodes and species in phylogenetic tree.
#' 
#' @references Gravel D., Canhan C.D., Beaudet M. and Messier C. Reconciling
#'     niche and neutrality: the continuum hypothesis. 2006. Ecology Letters 
#'     \doi{10.1111/j.1461-0248.2006.00884.x}
#'     
#' @export
#'
#' @examples
#'
calc_sbears <-
  function(x,
           phy,
           coords,
           method = c("single_site", "disp_assembly"),
           w_slope = 5,
           min_disp_prob=0.8,
           compute.node.by.sites = FALSE,
           make.node.label = FALSE
  ){
    # Enter and organize data:
    match <- picante::match.phylo.comm(phy, x)
    x <- match$comm
    
    # Extract species by nodes matrix with Herodotools
    if(make.node.label == TRUE){
      phy <- ape::makeNodeLabel(phy = phy, method = "number", prefix = "Node")
    }
    #spp_nodes <- t(get_spp_nodes(tree = phy, node.prefix = "Node")) # alternative
    
    # Run Ancestral Area Reconstruction:
    node.list <- list()
    node.anc.area.spat <- node.anc.area <- node.samp.mat <- matrix(NA, nrow = phy$Nnode, ncol = nrow(x), dimnames = list(phy$node.label, rownames(x)))
    
    for(i in 1:nrow(x)){
      node.list[[i]] <- phytools::fastAnc(phy, x[i, ])
      node.anc.area[, i] <- node.list[[i]]
      #print(i)
    }
    
    m.node.anc.area <- rowMeans(node.anc.area)
    sd.node.anc.area <- numeric()
    for(i in 1:nrow(node.anc.area)){
      sd.node.anc.area[i] <- sd(as.numeric(node.anc.area[i,]))
    }
    for(i in 1:nrow(node.anc.area)){
      # i = 1
      for(p in 1:ncol(node.anc.area)){
        # p = 1
        node.anc.area[i,p] <- pnorm(q = (node.anc.area[i, p] - m.node.anc.area[i])/sd.node.anc.area[i],
                                    mean=0,sd=1, lower.tail=TRUE)
      }
    }
    
    if (method == "disp_assembly"){
      r <- scales::rescale(geodist::geodist(x = coords, measure = "geodesic")/1000,diag = T, upper = T, c(0, 1)) # in km
      rownames(r)<-colnames(r)<-rownames(x)
      max_disp_dist<-sqrt(-log(min_disp_prob)/w_slope) # max distance
      anc_list<-list()
      site_values_pernode<-matrix(NA,nrow(x),nrow(x),dimnames=list(rownames(x),rownames(x)))
      for (i in 1:nrow(r)){
        r_below_threshold<-which(r[i,]<=max_disp_dist)
        r_pruned<-r[r_below_threshold,r_below_threshold]
        dist.decay<-matrix(NA, nrow = phy$Nnode, ncol = nrow(r_pruned), dimnames = list(phy$node.label, rownames(r_pruned)))
        
        for (k in rownames(r_pruned)){
          for (j in 1:nrow(node.anc.area)){
            for (p in colnames(r_pruned)){
              dist.decay[j,p]<-node.anc.area[j,k]*exp(1)^-(w_slope*r_pruned[k,p]^2)
            }
          }
        }
        anc_list[[i]]<-dist.decay
      }
      
      for (j in 1:nrow(node.anc.area)){
        for(l in 1:length(anc_list)){
          dist.decay<-anc_list[[l]]
          for (k in colnames(dist.decay)){
            site_values_pernode[l,k]<-anc_list[[l]][j,k]
          }
        }
        for(k in 1:ncol(node.anc.area)){
          site_values_pernode_site<-site_values_pernode[,k][which(!is.na(site_values_pernode[,k]==TRUE))]
          dens_site_node_pernode<-density(site_values_pernode_site,from=0,to=1)
          node.anc.area.spat[j,k]<-dens_site_node_pernode$x[which.max(dens_site_node_pernode$y)]
        }
      }
      node.anc.area<-node.anc.area.spat
    } else { node.anc.area=node.anc.area}
    
    # Compute a matrix of nodes by sites
    if(compute.node.by.sites==TRUE){
      comp_ada_nodes_sites <-
        function(phy, comm, long = FALSE){
          node_samp_mat <- matrix(NA,
                                  nrow = phy$Nnode,
                                  ncol = nrow(comm),
                                  dimnames = list(phy$node.label, rownames(comm)))
          names_node <- phy$node.label
          for(i in 1:length(names_node)){
            for(j in 1:nrow(comm)){
              comm_samp <- comm[, which(comm[j, ] == 1)]
              samp_nodes <- picante::prune.sample(samp = comm_samp, phylo = phy)$node.label
              node_samp_mat[i, j] <- ifelse(names_node[i] %in% samp_nodes, 1, 0)
            }
          }
          node_samp_mat <- t(node_samp_mat)
          if(long == TRUE){
            node_samp_mat <- phyloregion::dense2long(node_samp_mat)
            return(node_samp_mat)
          } else{
            return(node_samp_mat)
          }
        }
      
      node.samp.mat <- comp_ada_nodes_sites(phy = phy, comm = x, long = FALSE)
    } else {node.samp.mat<-"Nops..."}
    
    list_res <- list(reconstruction=node.anc.area, PD_nodes_by_sites = node.samp.mat)
    return(list_res)
  }