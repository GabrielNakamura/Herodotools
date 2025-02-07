
#' Computation of core ADA components
#'
#' @param x community matrix. Sites are rows and species are columns
#' @param phy phylogenetic tree as newick object
#' @param type character indicating the type of reconstruction to be used. Default is continuous
#' @param lik.threshold Logical, default is FALSE
#' @param threshold scalar indicating the threshold to be used when considering the presence of a lineage in a site
#' @param compute.node.by.sites Logical, TRUE (default) computes a matrix of node occurrence by site
#' @param make.node.label Logical, if TRUE (default) the nodes of the phylogeny will be named as the letter "N" preceding node number
#'
#' @return a list with three elements. reconstruction is the result of ancestral area reconstruction; phylogeny is the matrix containing
#'     the occurrence of nodes in sites and joint.phylo.obs is the joint occurrence of nodes and species in phylonetic tree
#' @export
#'
#' @examples
#' 
sbears <- 
  function(x, 
           phy, 
           type = c("discrete", "continuous"),
           lik.threshold = FALSE,
           threshold = 0.20, 
           compute.node.by.sites = TRUE, 
           make.node.label = TRUE
           ){
    # Enter and organize data:
    match <- picante::match.phylo.comm(phy, x)
    x <- match$comm
    
    root.age <- max(cophenetic(phy))/2
    
    # Extract species by nodes matrix with Herodotools
    if(make.node.label == TRUE){
      phy <- ape::makeNodeLabel(phy = phy, method = "number", prefix = "Node")
    }
    spp_nodes <- t(get_spp_nodes(tree = phy, node.prefix = "Node")) # alternative
    
    # Run Ancestral Area Reconstruction:
    node.list <- list()
    node.anc.area <- node.samp.mat <- matrix(NA, nrow = phy$Nnode, ncol = nrow(x), dimnames = list(phy$node.label, rownames(x)))
    
    if(type=="continuous"){
      for(i in 1:nrow(x)){
        node.list[[i]] <- phytools::fastAnc(phy, x[i, ])
        node.anc.area[, i] <- node.list[[i]]
        print(i)
      }
      if(lik.threshold == FALSE){
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
      } else {
        for(i in 1:nrow(x)){
          for (j in 1:phy$Nnode){
            ifelse(node.list[[i]][j] >= threshold, 1, 0)
          }
        }
        node.anc.area[,i]<-node.list[[i]]
      }
      
    } else {
      node.list.single<-list()
      for(i in 1:nrow(x)){
        node.list.single[[i]]<-apply(ape::ace(t(x)[,i], phy,type=type)$lik.anc,1,which.max)
        node.anc.area[,i]<-ifelse(node.list.single[[i]]==1,0,1) #
      }
    }
    
    # Compute a matrix of nodes by sites
    if(compute.node.by.sites==TRUE){
      node.samp.mat <- comp_ada_nodes_sites(phy = phy, comm = x, long = FALSE)
    } else {node.samp.mat<-"Nops..."}
    
    # Define joint occurrence matrix:
    joint.x <- as.matrix(cbind(x, t(node.anc.area))) # nodes and species occurrences in communities
    joint_x_occ <- ifelse(joint.x >= threshold, 1, 0)
    list_res <- list(reconstruction = node.anc.area, phylogeny = node.samp.mat, joint.phylo.obs = joint_x_occ)
    return(list_res)
  }
