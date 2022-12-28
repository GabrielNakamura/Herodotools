#' Get node ranges from BioGeoBEARS biome reconstruction model
#'
#' @description This is an auxiliary function that can be with a ancestral area reconstruction
#'    from BioGeoBEARS to obtain a data frame containing the ancestral area of occurrence (biome range)
#'    for all ancestors (nodes) 
#' @param BioGeoBEARS.data An object containing the result of BioGeoBEARS reconstruction
#' @param phyllip.file Phylip file, the same used in BioGeoBEARS package
#' @param tree A phylogenetic tree
#' @param max.range.size Scalar indicating the maximum number of biomes in which each node can occupy 
#'
#' @return Data frame with two columns, one indicating the node and other the ancestral biome range
#' 
#' @export
#'
#' @examples
#' data(resDEC) # ancestral area reconstruction from biogeobears
#' data(akodon_newick)
#' extdata_dir = system.file("extdata", package="Herodotools")
#' fn = paste(extdata_dir, "/geo_area_akodon.data", sep="")
#' node_area <- get_node_range_BioGeoBEARS(resDEC,phyllip.file = fn, akodon_newick,
#'                                         max.range.size = 4) # obtaining area for each node
#' 
get_node_range_BioGeoBEARS <-
  function(BioGeoBEARS.data,
           phyllip.file,
           tree,
           max.range.size){
    
    pkg_req <- c("BioGeoBEARS", "cladoRcpp")
    
    for(pkg in pkg_req) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
          paste0("Package '", pkg, "' must be installed to use this function."),
          call. = FALSE
        )
      }
    }
    
    tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = phyllip.file)
    trtable <- BioGeoBEARS::prt(tree, printflag=FALSE)
    prob_state <- BioGeoBEARS.data$ML_marginal_prob_each_state_at_branch_bottom_below_node
    areas <- BioGeoBEARS::getareas_from_tipranges_object(tipranges)
    states_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(areas = areas, 
                                                                    maxareas = max.range.size, 
                                                                    include_null_range = TRUE)
    ranges_list = NULL
    for (i in 1:length(states_list_0based)) {
      if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
      {
        tmprange = "_"
      } else {
        tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
      }
      ranges_list = c(ranges_list, tmprange)
    }
    
    range_probabilities <- as.data.frame(prob_state)
    row.names(range_probabilities) <- trtable$node
    names(range_probabilities) <- ranges_list
    teste_max <- apply(range_probabilities,
                       MARGIN = 1, 
                       function(i) colnames(range_probabilities)[which(i == max(i))
                       ]
    ) # taking the max probability value
    
    nodes.area <- data.frame(node = 1:length(unlist(teste_max)), biome = unlist(teste_max))
    nodes.range <- data.frame(area = nodes.area[-c(1:length(tree$tip.label)), 2])
    rownames(nodes.range) <- paste("N", 
                                    nodes.area[-c(1:length(tree$tip.label)), 1],
                                    sep = "")
    
    return(nodes.range)
    
  }