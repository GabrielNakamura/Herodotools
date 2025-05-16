plot_biogeobears_stochastic_map <- function(resDEC, bsm, map_index = 1, show.tip.label = TRUE, lwd = 5, label.offset = 0.1) {
  requireNamespace("BioGeoBEARS")
  requireNamespace("cladoRcpp")
  
  # Check map index
  n_maps <- length(bsm$RES_clado_events_tables)
  if (map_index < 1 || map_index > n_maps) {
    stop(paste0("Invalid 'map_index'. Must be between 1 and ", n_maps, "."))
  }
  
  # Setup
  scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
  stratified <- FALSE
  include_null_range <- TRUE
  
  if (is.null(resDEC$inputs$trfn) || is.null(resDEC$inputs$geogfn)) {
    stop("resDEC must contain 'inputs$trfn' and 'inputs$geogfn' paths to tree and geography files.")
  }
  
  max_range_size <- resDEC$inputs$max_range_size
  clado_events <- bsm$RES_clado_events_tables[[map_index]]
  ana_events <- bsm$RES_ana_events_tables[[map_index]]
  
  # Tip ranges and area states
  tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = resDEC$inputs$geogfn)
  areanames <- names(tipranges@df)
  
  states_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(
    areas = areanames,
    maxareas = max_range_size,
    include_null_range = include_null_range
  )
  
  colors_list <- get_colors_for_states_list_0based(
    areanames = areanames,
    states_list_0based = states_list_0based,
    max_range_size = max_range_size,
    plot_null_range = TRUE
  )
  
  # Modify results object with stochastic map
  resmod <- stochastic_map_states_into_res(
    res = resDEC,
    master_table_cladogenetic_events = clado_events,
    stratified = stratified
  )
  
  # Plot tree with text
  plot_BioGeoBEARS_results(
    results_object = resmod,
    analysis_titletxt = paste("Stochastic map", map_index),
    addl_params = list("j"),
    label.offset = label.offset,
    plotwhat = "text",
    cornercoords_loc = scriptdir,
    root.edge = TRUE,
    colors_list_for_states = colors_list,
    skiptree = FALSE,
    show.tip.label = show.tip.label
  )
  
  # Paint branches
  paint_stochastic_map_branches(
    res = resmod,
    master_table_cladogenetic_events = clado_events,
    colors_list_for_states = colors_list,
    lwd = lwd,
    lty = par("lty"),
    root.edge = TRUE,
    stratified = stratified
  )
  
  # Plot text again
  plot_BioGeoBEARS_results(
    results_object = resmod,
    analysis_titletxt = paste("Stochastic map", map_index),
    addl_params = list("j"),
    plotwhat = "text",
    cornercoords_loc = scriptdir,
    root.edge = TRUE,
    colors_list_for_states = colors_list,
    skiptree = TRUE,
    show.tip.label = show.tip.label
  )
  
  invisible(resmod)
}
