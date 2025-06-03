#' Get the node area states from the biogeographical stochastic mapping's result
#' 
#' @param bsm results from the function [Herodotools::calc_bsm()] 
#' @param BioGeoBEARS.data a BioGeoBEARS result model object 
#' @param phyllip.file path to the phyllip file used in the original model 
#' @param tree.path path for the phylogenetic tree used in the original model
#' @param max.range.size maximum range size used in the biogeographical reconstruction
#'
#' @returns a list with same length as the number of stochastic mapping. Each 
#'  list element is a single column data frame with node area states
#'
#' @export
#'
#' 
get_bsm_node_area <- function(
    bsm, 
    BioGeoBEARS.data,
    phyllip.file,
    tree.path,
    max.range.size 
){
  
  res <- BioGeoBEARS.data
  res$inputs$trfn <- tree.path
  res$inputs$geogfn <- phyllip.file
  
  clado_events_tables <- bsm$RES_clado_events_tables
  
  l_bsm_realizations <- purrr::map(1:length(clado_events_tables), function(i){
    
    # Convert the BSM into a modified res object
    master_table_cladogenetic_events = clado_events_tables[[i]]
    
    resmod = BioGeoBEARS::stochastic_map_states_into_res(
      res = res, 
      master_table_cladogenetic_events = master_table_cladogenetic_events
    )
    
    phy <- ape::read.tree(tree.path)
    
    node_area_BSM <- Herodotools::get_node_range_BioGeoBEARS(
      resmod,
      phyllip.file = phyllip.file,
      phy,
      max.range.size = max.range.size 
    )                                  
    
    node_area_BSM                              
  })
}