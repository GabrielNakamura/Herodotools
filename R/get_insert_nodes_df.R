get_insert_nodes_df <- function(bsm, phyllip.file, max.range.size) {
  ana_events <- bsm$RES_ana_events_tables[[1]]
  clado_events <- bsm$RES_clado_events_tables[[1]]
  
  # Step 1: Get area/state mappings using BioGeoBEARS internals
  tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = phyllip.file)
  areas <- BioGeoBEARS::getareas_from_tipranges_object(tipranges)
  states_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(
    areas = areas,
    maxareas = max.range.size,
    include_null_range = TRUE
  )
  
  ranges_list <- purrr::map_chr(states_list_0based, function(state) {
    if (length(state) == 1 && is.na(state)) {
      "_"
    } else {
      paste0(areas[state + 1], collapse = "")
    }
  })
  
  # Step 2: Extract post-split pseudo-events
  clado_shifts <- purrr::map_dfr(seq_len(nrow(clado_events)), function(i) {
    row <- clado_events[i, ]
    parent_state <- row$sampled_states_AT_nodes
    data <- list()
    
    for (desc in c("left_desc_nodes", "right_desc_nodes")) {
      desc_node <- row[[desc]]
      
      if (!is.null(desc_node) && !is.na(desc_node)) {
        
        # define post-split state
        if(desc == "left_desc_nodes") desc_state <- row$samp_LEFT_dcorner
        if(desc == "right_desc_nodes") desc_state <- row$samp_RIGHT_dcorner
        
        if (desc_state != parent_state) {
          data[[length(data) + 1]] <- data.frame(
            node = desc_node,
            ancestor = row$node,
            event_time = 1e-20,
            from_state = parent_state,
            to_state = desc_state,
            from_range = ranges_list[parent_state],
            to_range = ranges_list[desc_state],
            event_type = "post_split",
            source = "post_split"
          )
        }
      }
    }
    if (length(data) > 0) do.call(rbind, data) else NULL
  })
  
  # Step 3: True anagenetic events (already use range strings)
  ana_shifts <- ana_events %>%
    dplyr::filter(!is.na(event_type)) %>%
    dplyr::transmute(
      node = nodenum_at_top_of_branch,
      ancestor = ancestor,
      event_time = event_time,
      from_state = NA,
      to_state = NA,
      from_range = current_rangetxt,
      to_range = new_rangetxt,
      event_type = event_type,
      source = "anagenetic"
    )
  
  # Step 4: Combine both
  out <- dplyr::bind_rows(ana_shifts, clado_shifts)
  return(out)
}

