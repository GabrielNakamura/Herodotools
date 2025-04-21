get_bsm_node_area <- function(
    result_bsm = bsm_res, 
    BioGeoBEARS.data = resDEC,
    phyllip.file = phyllip_file,
    tree.path = "akodon.new",
    max.range.size = 4 
){
  
  res <- BioGeoBEARS.data
  res$inputs$trfn <- tree.path
  res$inputs$geogfn <- phyllip.file
  
  clado_events_tables <- result_bsm$RES_clado_events_tables
  
  l_bsm_realizations <- lapply(1:length(clado_events_tables), function(i){
    
    # Convert the BSM into a modified res object
    master_table_cladogenetic_events = clado_events_tables[[i]]
    
    resmod = stochastic_map_states_into_res(
      res=res, 
      master_table_cladogenetic_events=master_table_cladogenetic_events
    )
    
    phy <- read.tree(tree.path)
    
    node_area_BSM <- Herodotools::get_node_range_BioGeoBEARS(
      resmod,
      phyllip.file = phyllip.file,
      phy,
      max.range.size = max.range.size 
    )                                  
    
    node_area_BSM                              
  })
}