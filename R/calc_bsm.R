#' calculate Biogeographical Stochastic Mapping (BSM) with BioGeoBEARS
#' 
#' This is a wrapper function to the BioGeoBEARS::runBSM(). It uses the result of a BioGeoBEARS model 
#' as input to generate realized range area changes across the tree nodes and branches.
#' This mapping is done stochastically based on the model parameters. 
#'
#' @param BioGeoBEARS.data a BioGeoBEARS result model object 
#' @param phyllip.file path to the phyllip file used in the original model 
#' @param tree.path path for the phylogenetic tree used in the original model
#' @param max.maps maximum number of stochastic maps to try to generate 
#'        (passed to `maxnum_maps_to_try` in `runBSM()` from `{BioGeoBEARS}`)
#' @param n.maps.goal Target number of successfully generated stochastic maps
#' @param seed Optional seed for reproducibility of stochastic simulations.
#' @param save_after_every_try Logical. If `TRUE`, results will be saved after each attempt (default is `FALSE`).
#' @param ... Additional arguments passed to `runBSM()` from `{BioGeoBEARS}`
#'
#' @return A list (BSM output) containing the simulated mappings, including ancestral range estimates
#'         and event counts. This output is directly produced by `runBSM()` from `{BioGeoBEARS}`.
#'
#' @examples
#' \dontrun{
#' data("akodon_newick")
#' data("resDEC")  # BioGeoBEARS model result
#'
#' bsm_result <- calc_bsm(
#'   BioGeoBEARS.data = resDEC,
#'   phyllip.file = system.file("extdata", "geo_area_akodon.data", package = "Herodotools"),
#'   tree.path = system.file("extdata", "akodon.new", package = "Herodotools"),
#'   max.maps = 50,
#'   n.maps.goal = 10,
#'   seed = 1234
#' )
#' }
#'
#' @references
#' Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Frontiers of Biogeography*, 5(4).
#' 
#' Dupin, J., et al. (2017). Bayesian estimation of the global biogeographical history of the Solanaceae. *Journal of Biogeography*, 44(4), 887–899.
#'
#' @export


calc_bsm <- function(
    BioGeoBEARS.data,
    phyllip.file,
    tree.path,
    max.maps = 100, 
    n.maps.goal = 50,
    seed = NULL,
    save_after_every_try = FALSE,
    ...
) {
  
  
  res <- BioGeoBEARS.data
  res$inputs$trfn <- tree.path
  res$inputs$geogfn <- phyllip.file
  
  # Suppress cat() and message() from get_inputs_for_stochastic_mapping
  stochastic_mapping_inputs_list <- withCallingHandlers(
    {
      temp <- NULL
      invisible(capture.output({
        temp <- BioGeoBEARS::get_inputs_for_stochastic_mapping(res = res)
      }))
      temp
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  
  # Use temporary directory to avoid saving files to disk
  temp_savedir <- tempfile("bsm_tmp_")
  dir.create(temp_savedir)
  
  # Suppress cat() and message() from runBSM
  BSM_output <- withCallingHandlers(
    {
      temp <- NULL
      invisible(capture.output({
        temp <- BioGeoBEARS::runBSM(
          res,
          stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
          maxnum_maps_to_try = max.maps,
          nummaps_goal = n.maps.goal,
          seedval = seed,
          save_after_every_try = save_after_every_try,
          savedir = temp_savedir,
          ...
        )
      }))
      temp
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  
  return(BSM_output)
}
