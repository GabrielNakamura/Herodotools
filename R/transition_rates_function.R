#' Tip-based metrics of trait evolution
#' 
#' @details This function calculates three metrics that represent the macroevolutionary dynamic of traits for each species in a phylogenetic tree. They are 
#'     Transition rates, last transition time and stasis time, all originally proposed in Luza et al (2021). *Transition rates* indicate how many
#'     times the ancestral character has changed over time. *Stasis time* indicates the maximum branch length (time interval)
#'     which the current tip-character was maintained across the whole phylogeny. Finally, *last transition time* is the sum of branch lengths
#'     from the tip to the previous node with a reconstructed character equal to the current tip-character. Originally the function
#'     uses a stochastic character mapping on discrete trait data, but other types of trait reconstruction can be used.
#'     
#' @param tree A phylogenetic tree as an object of class "phylo"
#' @param trait A named vector containing the tip states for a discretely valued character. The names must match the tip labels of the tree.
#' @param nsim Number of simulations to stochastic character mapping.
#' @param method Tip-based metric, partial match to "transition_rates", "last_transition_time" and "stasis_time".
#' @return A list (length equal to nsim) with tip-based metrics estimated per species.
#' 
#' @author André Luza and Vanderlei Debastiani
#' 
#' @importFrom stats setNames
#' 
#' @example 
#' \dontrun{
#' data(rodent.phylo) # phylogenetic tree
#' data(trait) # categorical traits on species diet
#' trans_rates <- tip_based_trait_evo(tree=match_data$phy,trait =trait , # need to be named to worknsim = 50,method = c("transition_rates", "last_transition_time", "stasis_time"))
#' }
#' 
#' @export
#' 
tip_based_trait_evo <- function(tree,
                                trait, 
                                nsim = 1,
                                method = c("transition_rates", "last_transition_time", "stasis_time")
                                ) {
 
  pkg_req <- c("daee")
  
  for(pkg in pkg_req) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        paste0("Package '", pkg, "' must be installed to use this function."),
        call. = FALSE
      )
    }
  }
  
   ## Internal function
  adjacency.tree <- function(tree){
    temp <- ape::prop.part(tree)
    result <- matrix(0, nrow = length(tree$tip), ncol = length(temp), dimnames = list(tree$tip.label, tree$node.label))
    for(i in 1:ncol(result)){
      result[temp[[i]],i] <- 1
    }
    return(result)	
  }
  
  tree <- daee::node.tree(tree)$tree
  n.sp <- ape::Ntip(tree)
  n.no <- ape::Nnode(tree)
  spxnode <- adjacency.tree(tree)
  spxnode.edge <- matrix(NA, n.sp, n.no, dimnames = list(rownames(spxnode), colnames(spxnode)))
  spxnode.length <- matrix(NA, n.sp, n.no, dimnames = list(rownames(spxnode), colnames(spxnode)))
  for(i in 1:n.no){
    temp <- daee::tree.label.info(tree, colnames(spxnode)[i])
    spxnode.edge[spxnode[,i]==1, i] <- temp$edge
    spxnode.length[spxnode[,i]==1, i] <- temp$edge.length
  }
  sp.edge <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "Edge"))
  sp.length <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "Length"))
  for(i in 1:n.sp){
    temp <- daee::tree.label.info(tree, rownames(spxnode)[i])
    sp.edge[i,1] <- temp$edge  
    sp.length[i,1] <- temp$edge.length  
  }
  ## stochastic mapping of discrete traits via Bayesian inference 
  result.scm <-  phytools::make.simmap(tree, trait, model = "SYM", type = "discrete", nsim = nsim, message = FALSE, pi = "estimated", Q = "mcmc")
  ## Internal function 
  f.int <- function(scm, spxnode, sp.edge, trait, n.sp, n.no, method){
    states <- setNames(sapply(scm$maps, function(x) names(x)[1]), scm$edge[, 1])
    states <- states[as.character(ape::Ntip(scm) + 1:scm$Nnode)]
    sp.node.trait <- sweep(spxnode, 2, STATS = cbind(states), function(x, z) ifelse(x==1, z, NA))
    METHOD <- c("transition_rates", "last_transition_time", "stasis_time")
    method <- pmatch(method, METHOD)
    if (any(is.na(method))) {
      stop("\n Invalid method \n")
    }
    result <- data.frame(row.names = rownames(spxnode))
    if(any(method==1)){
      transitions <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "transition"))
      for(i in 1:n.sp){
        base.temp <- trait[rownames(spxnode)[i]]
        n.tra.temp <- 0
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            if(sp.node.trait[i,j]!=base.temp){
              n.tra.temp <- n.tra.temp+1
              base.temp <- sp.node.trait[i,j]
            }
          }
        }
        transitions[i,1] <- n.tra.temp
      }
      total.nodes <- cbind(apply(spxnode, 1, sum))
      prop.transitions <- transitions/total.nodes
      colnames(prop.transitions) <- "prop.transitions"
      result$transitions <- transitions
      result$total.nodes <- total.nodes
      result$prop.transitions <- prop.transitions
    }
    if(any(method==2)){
      time.transition <- matrix(NA, n.sp, 1, dimnames = list(rownames(spxnode), "last.transition.time"))
      for(i in 1:n.sp){
        base.temp <- trait[rownames(spxnode)[i]]
        time.tra.temp <- sp.length[i,1]
        go <- TRUE
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            if(sp.node.trait[i,j]!=base.temp){
              go <- FALSE
            } else{
              if(go){
                time.tra.temp <- time.tra.temp+spxnode.length[i,j]
              }
            }
          }
        }
        time.transition[i,1] <- time.tra.temp
      }
      result$last.transition.time <- time.transition
    }
    if(any(method==3)){
      scm.maps.max <- sapply(scm$maps, function(x) x[which(x==max(x))], simplify = FALSE)
      stasis.time <- matrix(NA, n.sp, 1, dimnames=list(rownames(spxnode), "stasis.time"))
      for(i in 1:n.sp){
        time.sta.temp <- scm.maps.max[[sp.edge[i,1]]]
        for(j in n.no:1){
          if(spxnode[i,j]!=0){
            time.sta.temp <- max(time.sta.temp, scm.maps.max[[spxnode.edge[i,j]]])
          }
        }
        stasis.time[i,1] <- time.sta.temp
      }
      result$stasis.time <- stasis.time
    }
    return(result)
  }
  if(nsim==1){
    RES <- sapply(list(result.scm), f.int, spxnode = spxnode, sp.edge = sp.edge, trait = trait, n.sp = n.sp, n.no = n.no, method = method, simplify = FALSE)
  } else{
    RES <- sapply(result.scm, f.int, spxnode = spxnode, sp.edge = sp.edge, trait = trait, n.sp = n.sp, n.no = n.no, method = method, simplify = FALSE)
  }
  return(RES)
}
