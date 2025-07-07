calc_equal_splits <- function(tree, node_area, current_area = NULL, ignore_nodes = NULL) {
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  if (!is.data.frame(node_area)) stop("node_area must be a data.frame with rownames as node labels and a column 'area'")
  
  edge <- tree$edge
  edge.length <- tree$edge.length
  n_tips <- length(tree$tip.label)
  
  # Map node labels to node numbers
  node_labels <- setNames(n_tips + seq_along(tree$node.label), tree$node.label)
  
  # Determine nodes to ignore (e.g., inserted nodes like "new_N...")
  if (is.null(ignore_nodes)) {
    ignore_nodes <- grep("^new_N", tree$node.label, value = TRUE)
  }
  ignored_node_nums <- unname(node_labels[ignore_nodes])
  
  # Convert node_area to map from node number to area string
  node_area$node <- node_labels[rownames(node_area)]
  area_by_node <- setNames(as.character(node_area$area), node_area$node)
  
  # Initialize output ED vector
  ed <- setNames(numeric(n_tips), tree$tip.label)
  
  # Loop through each tip
  for (tip in 1:n_tips) {
    node <- tip
    decay_index <- 0
    pending_length <- 0
    ed_value <- 0
    
    while (node != (n_tips + 1)) {  # walk toward root
      parent_edge_index <- which(edge[, 2] == node)
      parent_node <- edge[parent_edge_index, 1]
    
      
      # If we are stopping due to current_area
      if (!is.null(current_area)) {
        parent_area <- area_by_node[as.character(parent_node)]
        if (is.na(parent_area) || !grepl(paste0("[", current_area, "]"), parent_area)) {
          # At stop point: include ignored node only now
          ed_value <- ed_value + (pending_length * (0.5 ^ decay_index))
          break
        }
      }
      
      branch_len <- edge.length[parent_edge_index]
      pending_length <- pending_length + branch_len
      
      # Normal case: accumulate ED only if not an ignored node
      if (!(parent_node %in% ignored_node_nums)) {
        ed_value <- ed_value + (pending_length * (0.5 ^ decay_index))
        pending_length <- 0
        decay_index <- decay_index + 1
      }
      
      node <- parent_node
    }
    
    ed[tip] <- ed_value
  }
  
  return(ed)
}


