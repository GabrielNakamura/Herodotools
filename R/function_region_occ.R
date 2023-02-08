#' Species occurrence in each evoregion
#'
#' Auxiliary function to produce an occurrence data frame of species in each region. This object can be 
#'     with \code{\link{get_tipranges_to_BioGeoBEARS}} function to produce a phyllip file needed in BioGeoBEARS
#'     ancestral area reconstruction
#'
#' @param comm Occurrence data frame with species in colums and rows corresponding to assemblages
#' @param site.region a character vector indicating the membership of each assemblage to regions
#' @importFrom rlang .data
#' @return An occurrence data frame with species as rownames and regions as columns
#' @export
#' 
#' @seealso \code{\link{get_tipranges_to_BioGeoBEARS}}
#'
#' @author Gabriel Nakamura and Arthur Rodrigues
#'
#' @examples
#' data(akodon_sites)
#' data(akodon_newick)
#' akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
#' spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
#' akodon_pa_tree <- akodon_pa[, spp_in_tree]
#' data(regions)
#' site_evoregion <- regions$Cluster_Evoregions
#' get_region_occ(comm=akodon_pa_tree,site.region=site_evoregion)
get_region_occ <- 
  function(comm, site.region){
  evoregion.data <- 
    # bind compostion matrix and the site evoregion
    dplyr::bind_cols(
      comm,
      site.region =  site.region
    ) %>% 
    # transform the matrix in a longer dataframe 
    tidyr::pivot_longer(
      cols = 1:ncol(comm), 
      names_to = "species",
      values_to = "presence"
    ) %>% 
    # filter only the rows with the presence of a species
    dplyr::filter(.data$presence == 1) %>% 
    # Count the number of sites in a region that a species is present
    dplyr::group_by(.data$species, site.region) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    # Calculate the propotional area of a species in each region
    dplyr::group_by(.data$species) %>% 
    dplyr::mutate(
      species.total = sum(.data$n), 
      prec.occupation = .data$n/.data$species.total
    ) %>% 
    dplyr::ungroup() %>% 
    # Filter for each species the evoregions with more than 25% of area ocupied 
    # by the species
    dplyr::filter(
      .data$prec.occupation >= .25
    ) %>% 
    dplyr::mutate(
      area = LETTERS[site.region]
    ) %>% 
    dplyr::select(.data$species, .data$area) %>% 
    dplyr::mutate(value = 1) %>% 
    # transform in a dataframe with species in rows and areas ocupied [0-1] in columns
    tidyr::pivot_wider(
      id_cols = .data$species,
      names_from = .data$area, 
      names_sort = T,
      values_from = .data$value,
      values_fill = 0
    ) %>% 
    as.data.frame()
  species.names <- evoregion.data[,1]
  a.regions <- evoregion.data[,-1]
  rownames(a.regions) <- species.names  
  return(a.regions)
}
