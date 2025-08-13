calc_node_density <- function(tree, ancestral.area, current.area = NULL) {
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  if (!is.data.frame(ancestral.area)) stop("ancestral.area must be a data.frame with rownames as node labels and a column 'area'")
  
  edge <- tree$edge
  edge.length <- tree$edge.length
  n_tips <- length(tree$tip.label)
  
  # Map node labels to node numbers
  node_labels <- setNames(n_tips + seq_along(tree$node.label), tree$node.label)
  
  # Identify inserted nodes to ignore
  ignore_nodes <- grep("^new_N", tree$node.label, value = TRUE)
  ignored_node_nums <- unname(node_labels[ignore_nodes])
  
  # Convert ancestral.area to map from node number to area
  ancestral.area$node <- node_labels[rownames(ancestral.area)]
  area_by_node <- setNames(as.character(ancestral.area$area), ancestral.area$node)
  
  # Compute node depths from root
  node_depths <- ape::node.depth.edgelength(tree)
  age_from_present <- abs(node_depths - node_depths[1])
  
  nd <- setNames(numeric(n_tips), tree$tip.label)
  
  for (tip in 1:n_tips) {
    node <- tip
    node_count <- 0
    stopping_node <- NULL
    
    while (node != (n_tips + 1)) {
      parent_edge_index <- which(edge[, 2] == node)
      parent_node <- edge[parent_edge_index, 1]
      
      # Check stopping condition
      if (!is.null(current.area)) {
        parent_area <- area_by_node[as.character(parent_node)]
        if (is.na(parent_area) || !grepl(paste0("[", current.area, "]"), parent_area)) {
          stopping_node <- node
          break
        }
      }
      
      if (!(parent_node %in% ignored_node_nums)) {
        node_count <- node_count + 1
      }
      
      node <- parent_node
    }
    
    # Determine denominator: depth of stopping node or root
    if (is.null(stopping_node)) {
      stopping_node <- n_tips + 1  # root
    }
    depth <- age_from_present[stopping_node]
    
    nd[tip] <- node_count / depth
  }
  
  nd[is.na(nd)] <- 0
  
  return(nd)
}
