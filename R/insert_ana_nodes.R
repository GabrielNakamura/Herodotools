insert_nodes <- function(tree, inserts, node_area = NULL) {
  library(ggtree)
  library(dplyr)
  library(tibble)
  
  # Handle multiple insert tables via recursion
  if (is.list(inserts) && all(purrr::map_lgl(inserts, is.data.frame))) {
    if (!is.null(node_area)) {
      if (!is.list(node_area) || length(node_area) != length(inserts)) {
        stop("When 'inserts' is a list, 'node_area' must also be a list of the same length.")
      }
    } else {
      node_area <- vector("list", length(inserts))
    }
    
    results <- purrr::map2(inserts, node_area, ~insert_nodes(tree, .x, .y))
    return(results)
  }
  
  # Ensure node labels are present on the tree
  if (is.null(tree$node.label)) {
    tree$node.label <- paste0("N", seq_len(tree$Nnode) + Ntip(tree))
  }
  
  # Build initial ggtree data
  p <- ggtree(tree)
  gdata <- p$data
  
  # Ensure child and parent exist in the tree
  if (!all(inserts$child %in% gdata$node)) {
    stop("Some 'child' values in inserts do not exist in the tree.")
  }
  if (!all(inserts$parent %in% gdata$node)) {
    stop("Some 'parent' values in inserts do not exist in the tree.")
  }
  
  # Check that node labels match rownames in node_area
  if (!is.null(node_area)) {
    expected_nodes <- tree$Nnode
    actual_nodes <- nrow(node_area)
    
    if (expected_nodes != actual_nodes) {
      stop(paste0("Node count mismatch: tree has ", expected_nodes,
                  " internal nodes, but node_area has ", actual_nodes, " rows."))
    }
    
    if (!all(tree$node.label == rownames(node_area))) {
      stop("Mismatch between tree node labels and rownames in node_area.")
    }
  }
  
  # Ensure child-parent relationship is valid in the tree
  valid_edges <- gdata %>% select(parent, node) %>% filter(parent != 0) %>% suppressMessages()
  insert_edges <- inserts %>% select(parent, child)
  if (!all(paste(insert_edges$parent, insert_edges$child) %in%
           paste(valid_edges$parent, valid_edges$node))) {
    stop("Some insertions refer to invalid parent-child relationships.")
  }
  
  if (any(duplicated(inserts %>% select(parent, child, event_time)))) {
    stop("Duplicated insertion times for the same branch are not allowed.")
  }
  
  new_data <- gdata
  
  # Add node_area info to ggtree data (for original nodes)
  if (!is.null(node_area)) {
    node_lookup <- rownames(node_area)
    node_area_values <- node_area$area
    label_to_area <- setNames(node_area_values, node_lookup)
    new_data$node_area <- label_to_area[paste0("N", new_data$node)]
  }
  
  max_node <- max(new_data$node)
  inserts <- inserts %>% arrange(child, event_time)
  
  for (child_id in unique(inserts$child)) {
    child_inserts <- inserts %>% filter(child == child_id) 
    original_child_row <- new_data %>% filter(node == child_id) %>% suppressMessages()
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
    
    last_event_time <- tail(child_inserts$event_time, 1)
    final_bl <- original_bl - last_event_time
    
    new_data <- new_data %>%
      mutate(
        parent = ifelse(node == child_id, previous_node, parent),
        branch.length = ifelse(node == child_id, final_bl, branch.length)
      ) %>% 
      suppressMessages()
  }
  
  new_tree <- as.phylo(new_data) 
  
  new_node_area <- new_data %>% 
    as.data.frame() %>% 
    transmute(node = label, area = node_area) %>%
    filter(!is.na(area)) %>% 
    column_to_rownames("node")
  
  return(list(
    phylo = new_tree,
    node_area = new_node_area
  ))
}
