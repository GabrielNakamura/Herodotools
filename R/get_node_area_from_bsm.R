get_node_area_from_bsm <-function(bsm, res_model, include_null_range = TRUE) {
  if (!"RES_clado_events_tables" %in% names(bsm)) {
    stop("Input does not contain RES_clado_events_tables")
  }
  
  clado <- bsm$RES_clado_events_tables[[1]] %>% 
    filter(node.type != "tip")
  if (!"sampled_states_AT_nodes" %in% colnames(clado) || !"node" %in% colnames(clado)) {
    stop("Clado events table must contain 'node' and 'sampled_states_AT_nodes' columns")
  }
  
  # Get area names and max range size from original results
  tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(res_model$inputs$geogfn)
  areas <- BioGeoBEARS::getareas_from_tipranges_object(tipranges)
  max_range_size <- res_model$inputs$max_range_size
  
  # Get state-to-area mapping (0-based indices)
  states_list <- cladoRcpp::rcpp_areas_list_to_states_list(
    areas = areas,
    maxareas = max_range_size,
    include_null_range = include_null_range
  )
  
  # Convert state indices to range codes (like "AB", "C", etc.)
  state_to_range <- purrr::map_chr(states_list, function(indices) {
    if (length(indices) == 1 && is.na(indices)) {
      "_"
    } else {
      paste(areas[indices + 1], collapse = "")
    }
  })
  
  # Map numeric states to range codes
  node_labels <- paste0("N", clado$node)
  state_indices <- clado$sampled_states_AT_nodes   # convert from 0-based to 1-based
  valid_rows <- !is.na(state_indices) & state_indices <= length(state_to_range)
  
  node_area <- data.frame(
    area = state_to_range[state_indices[valid_rows]],
    row.names = node_labels[valid_rows],
    stringsAsFactors = FALSE
  )
  
  return(node_area)
}


