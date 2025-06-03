#' Get data frame of transition states from Biogeographical Stochastic Mapping
#'
#' @param bsm results from the function [Herodotools::calc_bsm()]
#' @param phyllip.file path to the phyllip file used in the original model
#' @param max.range.size maximum range size used in the biogeographical reconstruction
#' @param include_null_range default TRUE. 
#'
#' @returns a list with same length as the number of stochastic mapping. Each 
#'  list element is a data frame with columns:
#'    * `child` = child node number,
#'    * `parent` = ancestor node number,
#'    * `event_time` = time of the event from parent to child, 
#'    * `node_area` = the range area in the event,
#'    * `event_type` = the type of event,
#'    * `source` = can be 'cladogenesis' or 'anagenesis'
#'    
#' @export
#'

get_insert_df <- function(bsm, phyllip.file, max.range.size, include_null_range = TRUE) {
  n_maps <- length(bsm$RES_clado_events_tables)
  
  # Step 1: Get area/state mappings
  tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = phyllip.file)
  areas <- BioGeoBEARS::getareas_from_tipranges_object(tipranges)
  states_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(
    areas = areas,
    maxareas = max.range.size,
    include_null_range = include_null_range
  )
  ranges_list <- purrr::map_chr(states_list_0based, function(state) {
    if (length(state) == 1 && is.na(state)) {
      "_"
    } else {
      paste0(areas[state + 1], collapse = "")
    }
  })
  
  # Step 2: Loop over each stochastic map
  insert_dfs <- purrr::map(seq_len(n_maps), function(i) {
    ana_events <- bsm$RES_ana_events_tables[[i]]
    clado_events <- bsm$RES_clado_events_tables[[i]]
    
    # Post-split events
    clado_shifts <- purrr::map_dfr(seq_len(nrow(clado_events)), function(j) {
      row <- clado_events[j, ]
      parent_state <- row$sampled_states_AT_nodes
      data <- list()
      
      for (desc in c("left_desc_nodes", "right_desc_nodes")) {
        desc_node <- row[[desc]]
        if (!is.null(desc_node) && !is.na(desc_node)) {
          desc_state <- if (desc == "left_desc_nodes") row$samp_LEFT_dcorner else row$samp_RIGHT_dcorner
          if (desc_state != parent_state) {
            data[[length(data) + 1]] <- data.frame(
              child = desc_node,
              parent = row$node,
              event_time = 1e-20,
              node_area = ranges_list[desc_state],
              event_type = row$clado_event_type,
              source = "cladogenesis"
            )
          }
        }
      }
      if (length(data) > 0) do.call(rbind, data) else NULL
    })
    
    # Anagenetic events
    ana_shifts <- ana_events %>%
      dplyr::filter(!is.na(event_type)) %>%
      dplyr::transmute(
        child = nodenum_at_top_of_branch,
        parent = ancestor,
        event_time = event_time,
        node_area = new_rangetxt,
        event_type = event_type,
        source = "anagenesis"
      )
    
    dplyr::bind_rows(ana_shifts, clado_shifts)
  })
  
    return(insert_dfs)
}
