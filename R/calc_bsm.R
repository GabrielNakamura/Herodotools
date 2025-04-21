#' calculate Biogeographical Stochastic Mapping (BSM) with BioGeoBEARS
#'
#' @param BioGeoBEARS.data a BioGeoBEARS result model object 
#' @param phyllip.file path to the phyllip file 
#' @param tree.path path for the phylognetic tree
#' @param max.maps max number of stocastic maps
#' @param n.maps.goal number of desired stochastic maps
#' @param seed seed for reproducibility
#' @param save_after_every_try save BSM results on disk
#' @param ... other parameters passed to BioGeoBEARS::runBSM
#'
#' @returns
#' @export
#'
#' @examples
calc_bsm <- function(
    BioGeoBEARS.data,
    phyllip.file,
    tree.path,
    max.maps = 100, 
    n.maps.goal = 50,
    seed = NULL,
    save_after_every_try = FALSE,
    ...
){
  
  res <- BioGeoBEARS.data
  res$inputs$trfn <- tree.path
  res$inputs$geogfn <- phyllip.file
  
  stochastic_mapping_inputs_list <- BioGeoBEARS::get_inputs_for_stochastic_mapping(
    res = res
  )
  
  BSM_output <- BioGeoBEARS::runBSM(
    res,
    stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
    maxnum_maps_to_try = max.maps,
    nummaps_goal = n.maps.goal,
    seedval = seed,
    save_after_every_try = save_after_every_try,
    ...
  )
  
  BSM_output
  
}