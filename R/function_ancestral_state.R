#' Ancestral State per Assemblage
#'
#' @param tree A newick phylogenetic tree object
#' @param ancestral.area  A two column data frame. Lines are nodes the column are the biomes/region of occurrence
#'     for each ancestors. Can be obtained by using (\code{\link{get.node.range_BioGeoBEARS}}) 
#' @param prefix A single character string to be used to name nodes
#'
#' @return A data frame with assemblages in lines and nodes in columns. 
#'     Each cell contains the ancestral area/Ecoregion of occurrence for each node and its respective species. 
#'  
#'     
#' @author Gabriel Nakamura <gabriel.nakamura.souza@@gmail.com>
#'    
#' @export
#'
#' @examples
#' 
ancestral_state <- function(tree, 
                            ancestral.area, 
                            prefix = "N")
  {
  spxnode <- spp_nodes(tree = tree, node.prefix = prefix)
  AS <- matrix(data = 0, nrow = nrow(spxnode), 
               ncol = ncol(spxnode), 
               dimnames = list(rownames(spxnode), 
                               colnames(spxnode)
               )
  )
  for(i in 1:nrow(AS)){
    pres_nodes <- which(spxnode[i, ] == 1)
    absence_nodes <- which(spxnode[i, ] == 0)
    AS[i, pres_nodes] <- ancestral.area[i, ]
    AS[i, absence_nodes] <- NA
  }
  return(AS)
}
