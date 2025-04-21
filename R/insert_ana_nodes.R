insert_ana_nodes <- function(tree, inserts, node_area = NULL) {
  library(ggtree)
  library(dplyr)
  
  p <- ggtree(tree)
  gdata <- p$data
  new_data <- gdata 
  
  # add original node area info
  if(!is.null(node_area))  new_data$node_area <- node_area$area
  
  max_node <- max(gdata$node)
  
  # Sort inserts by child and increasing event_time (from parent toward tip)
  inserts <- inserts %>%
    arrange(child, event_time)
  
  for (child_id in unique(inserts$child)) {
    child_inserts <- inserts %>% filter(child == child_id)
    
    # Get original child row (to be reconnected at the end)
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
      new_row$x <- new_x
      new_row$branch.length <- event_time - (if (i == 1) 0 else child_inserts$event_time[i - 1])
      new_row$parent <- previous_node
      new_row$isTip <- FALSE
      new_row$node_area <- area_val
      
      new_data <- bind_rows(new_data, new_row)
      
      previous_node <- new_node
    }
    
    # Modify the original child to descend from the last inserted node
    last_event_time <- tail(child_inserts$event_time, 1)
    final_bl <- original_bl - last_event_time
    
    new_data <- new_data %>%
      mutate(
        parent = ifelse(node == child_id, previous_node, parent),
        branch.length = ifelse(node == child_id, final_bl, branch.length)
      )
  }
  
  # new node_area
  new_node_area <- new_data %>% 
    mutate(node = paste0("N", node)) %>% 
    select(node, node_area) %>% 
    rename(area = node_area) %>% 
      tibble::column_to_rownames("node")
    
  # Convert back to tree
  new_tree <- as.phylo(new_data)
  
  return(list(data = new_data, phylo = new_tree, node_area = new_node_area))
}
