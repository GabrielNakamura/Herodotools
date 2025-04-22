insert_ana_nodes <- function(tree, inserts, node_area = NULL) {
  library(ggtree)
  library(dplyr)
  library(tibble)
  
  # Assuming you have a phylo object `tree` and a node_area table
  expected_nodes <- tree$Nnode
  actual_nodes <- nrow(node_area)
  
  if (expected_nodes != actual_nodes) {
    stop(paste0("Node count mismatch: tree has ", expected_nodes,
                " internal nodes, but node_area has ", actual_nodes, " rows."))
  }
  
  # Ensure node labels are present on the tree
  if (is.null(tree$node.label)) {
    tree$node.label <- paste0("N", seq_len(tree$Nnode) + Ntip(tree))
  }
  
  # Check that node labels match rownames in node_area
  if (!is.null(node_area)) {
    expected_labels <- tree$node.label
    if (!all(tree$node.label == rownames(node_area))) {
      stop("Mismatch between tree node labels and rownames in node_area.")
    }
  }
  
  # Build initial ggtree data
  p <- ggtree(tree)
  gdata <- p$data
  new_data <- gdata
  
  # Add node_area info to ggtree data (for original nodes)
  if (!is.null(node_area)) {
    node_lookup <- rownames(node_area)
    node_area_values <- node_area$area
    label_to_area <- setNames(node_area_values, node_lookup)
    new_data$node_area <- label_to_area[paste0("N", new_data$node)]
  }
  
  max_node <- max(new_data$node)
  
  # Sort inserts to process each child in time order
  inserts <- inserts %>% arrange(child, event_time)
  
  for (child_id in unique(inserts$child)) {
    child_inserts <- inserts %>% filter(child == child_id)
    
    # Original child data
    original_child_row <- new_data %>% filter(node == child_id)
    original_parent <- original_child_row$parent
    original_bl <- original_child_row$branch.length
    original_x <- original_child_row$x
    
    previous_node <- original_parent
    
    for (i in seq_len(nrow(child_inserts))) {
      ins <- child_inserts[i, ]
      event_time <- ins$event_time
      area_val <- ins$node_area
      
      if (event_time >= original_bl || event_time <= 0) {
        stop("event_time must be > 0 and < original branch length")
      }
      
      # Create new node
      max_node <- max_node + 1
      new_node <- max_node
      new_x <- original_x - (original_bl - event_time)
      
      new_row <- original_child_row
      new_row$node <- new_node
      new_row$label <- paste0("ana_N", new_node)
      new_row$x <- new_x
      new_row$branch.length <- event_time - (if (i == 1) 0 else child_inserts$event_time[i - 1])
      new_row$parent <- previous_node
      new_row$isTip <- FALSE
      new_row$node_area <- area_val
      
      new_data <- bind_rows(new_data, new_row)
      previous_node <- new_node
    }
    
    # Adjust original child to attach to last inserted node
    last_event_time <- tail(child_inserts$event_time, 1)
    final_bl <- original_bl - last_event_time
    
    new_data <- new_data %>%
      mutate(
        parent = ifelse(node == child_id, previous_node, parent),
        branch.length = ifelse(node == child_id, final_bl, branch.length)
      )
  }
  
  # Convert to phylo
  new_tree <- as.phylo(new_data)
  
  # Reorder node_area using helper
  new_node_area <- match_node_area(new_tree, new_data %>% 
                                     transmute(node = paste0("ana_N", node), area = node_area) %>%
                                     filter(!is.na(area)))
  
  return(list(
    phylo = new_tree,
    node_area = new_node_area
  ))
}



match_node_area <- function(phylo, node_area) {
  internal_nodes <- phylo$node.label
  ana_nodes <- grep("^ana_N", rownames(node_area), value = TRUE)
  ordered_labels <- c(phylo$tip.label, internal_nodes, ana_nodes)
  
  matched <- node_area[ordered_labels[ordered_labels %in% rownames(node_area)], , drop = FALSE]
  return(matched)
}

