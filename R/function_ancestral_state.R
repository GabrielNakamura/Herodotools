
#' Ancestral State per Assemblage
#'
#' @param tree A newick phylogenetic tree object
#' @param ancestral.area  A two column data frame. Lines are nodes the columns are the biomes/region of occurrence
#'     for each ancestors. Can be obtained by using (\code{\link{get.node.range_BioGeoBEARS}}) 
#' @param prefix A single character string to be used to name nodes
#'
#' @return A data frame with assemblages in lines and nodes in columns. 
#'     Each cell contains the ancestral area of occurrence for each node. 
#'     
#' @export
#'
#' @examples
ancestral_state <- function(tree, 
                            ancestral.area, 
                            prefix = "N")
  {
  spxnodes <- spp_nodes(tree = tree, node.prefix = prefix)
  AS <- matrix(data = 0, nrow = nrow(spxnode), 
               ncol = ncol(spxnode), 
               dimnames = list(rownames(spxnode), 
                               colnames(spxnode)
               )
  )
  spxnodes_area <- cbind(spxnodes, ancestral.area)
  AS_res <- apply(spxnodes_area, MARGIN = 1, function(x){
    AS[which(x == 1), ] <- as.character(x[ncol(spxnodes_area)])
    AS[which(x == 0), ] <- NA
  })

  return(AS_res)
  
}
