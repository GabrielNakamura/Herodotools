#' Proportional contribution of an ancestral area to the species in other region
#' 
#' 
#' For the species present in a site, the function calculates the proportional 
#'  contribution of ancestral areas to dispersal of lineage in the site's region. 
#'
#' @param W Occurrence matrix, rows are assemblages and columns are species 
#' @param tree Phylogenetic tree in newick format
#' @param ancestral.area One column data frame with nodes in rows and one column indicating the occurrence (biome/ecoregion) area of nodes
#' @param biogeo One column data frame with assemblages in rows and their respective biome/ecoregion
#'
#' @author Arthur V Rodrigues <rodrigues.arthur.v@@gmail.com>
#'
#' @return A data frame with assemblages in the rows and regions in the columns. The values indicates
#'     the percentage of contribution of each region to each assemblage. NA represents no contribution
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' data(akodon_sites) # occurrence matrix
#' akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
#' data(akodon_newick) # phylogenetic tree
#' spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
#' akodon_pa_tree <- akodon_pa[, spp_in_tree]
#' data(regions) # biogeographic region
#' data(resDEC) # output from ancestral area reconstruction
#' node.area <- get_node_range_BioGeoBEARS(resDEC,
#'                                         phyllip.file = here::here("inst", 
#'                                         "extdata", "geo_area_akodon.data")
#'                                         ,akodon_newick,max.range.size = 4) 
#' calc_dispersal_from(W=akodon_pa_tree,
#'                     tree=akodon_newick,
#'                     ancestral.area=node.area,biogeo=regions) # historical dispersal analysis
#' }
#' 
#' @seealso \code{\link{get_node_range_BioGeoBEARS}}
#' 
calc_dispersal_from <- function(W,
                           tree,
                           ancestral.area, 
                           biogeo){
  
  nodes.list <- get_nodes_info_core(W = W, tree = tree, ancestral.area = ancestral.area, biogeo = biogeo)
  AS <- calc_ancestral_state(tree = tree, ancestral.area = ancestral.area, prefix = "N")
  names_spComm <- colnames(W)
  
  #matrix to receive the area from which dispersal ocurred
  dispersal_from <-  matrix("-",
                            nrow = nrow(W),
                            ncol = ncol(W),
                            dimnames = list(rownames(W), colnames(W)))
  
  ## NA means the dispersal was previous to the root node of phylogenetic tree used 
  for(site in 1:length(nodes.list)){
    for(sp in 1:length(nodes.list[[site]]$disp.anc.node)){
      pres <- which(W[site, ] >= 1)
      pres <- names_spComm[pres]
      node <- nodes.list[[site]]$disp.anc.node[sp]
      node.name  <- paste0("N", node)
      if(node.name == "NNA"){ dispersal_from[site, pres[sp]] <- NA
      }else{
        dispersal_from[site, pres[sp]] <- AS[node.name, pres[sp]]}
    }
  }
  
  l.freq.area <- lapply(1:nrow(dispersal_from), function(i){
    pres <- dispersal_from[i, ][W[i, ] >= 1]
    table(pres)/sum(table(pres))
  })
  
  areas <- unique(ancestral.area[ ,1])
  disp_from_frequency <- matrix(NA, nrow = nrow(W), ncol = length(areas))
  rownames(disp_from_frequency) <- row.names(W)
  colnames(disp_from_frequency) <- areas
  
  for(i in 1:nrow(disp_from_frequency)){
    temp.freq<- l.freq.area[[i]]
    temp.col <- colnames(disp_from_frequency) %in% names(temp.freq)
    
    disp_from_frequency[i, temp.col] <- temp.freq
  }
  
  return(disp_from_frequency)
  
}

