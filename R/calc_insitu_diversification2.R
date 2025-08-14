calc_insitu_diversification2 <- function(W,
                                         tree,
                                         ancestral.area, 
                                         biogeo, 
                                         diversification = "jetz",
                                         type = "equal.splits"){
  
  
  if(all(diversification != c("jetz", "freckleton")) == TRUE){
    stop("Tip-based diversification measures must be one of 'jetz' or 'freckleton'")
  }
  
  # Function to compute harmonic mean including zeros
  # x is a vector of a site
  # NA are removed from the calculation 
  harmonic_mean <- function(x) {
    x <- x[!is.na(x)]
    if (any(x == 0)) return(0)  # harmonic mean is 0 if any value is 0
    length(x) / sum(1 / x)
  }
  
 
  ## step 1 - calc tip rates for all species in all areas
  
  # define the unique areas 
  unique_ranges_colapsed <- paste(unique(biogeo[,1]), collapse = "")
  unique_areas <- unique(unlist(strsplit(unique_ranges_colapsed, "")))
  names(unique_areas) <- unique_areas
  
  
  if (diversification == "jetz" ){
    
    ed_total <- calc_ed(tree, ancestral.area, type = type)
    
    l_ed_partial <- lapply(unique_areas, function(area){
      ed_local <- calc_ed(tree, ancestral.area, current.area = area, type = type)
      
      ed_local/ed_total
    })
    
    jetz_total <- 1/ed_total
    l_jetz_partial <- lapply(l_ed_partial, function(ed_area){
      jetz_total*ed_area
    })
    
    # Initialize an empty matrix with the same dimensions and names as W
    matrix_XJetz <- W
    
    # Loop through each community
    for (i in seq_len(nrow(W))) {
      region <- biogeo[i, 1]
      jetz_vals <- l_jetz_partial[[region]][colnames(W)]
      pres_only <- ifelse(W[i, ] == 0, NA, W[i, ])
      matrix_XJetz[i, ] <- pres_only * jetz_vals
    }
    
    
    
    matrix_XJetz <- ifelse(matrix_XJetz == 0, 0.00001, matrix_XJetz)

    
    # Apply to each row
    JetzLocalComm_harmonic <- apply(matrix_XJetz, 1, harmonic_mean)
    
    
    
 
    
    
    
  }
  
  # compute total diversification 
  #look at type and compute the ed metric for 
  # using each area as the current area per iteration. 
  
  # if freckleton:
  #  # compute total diversification 
  # compute the nd for all the areas as current area 
  
  ## step 2 - compute the mean rates per site
  
  ## calculating harmonic mean for Jetz
  
  # sum of communities local diversification
  
  #sum of communities total diversification
  #harmonic mean for Jetz total diversification
  #harmonic mean for Jetz local diversification
  
  
  #substitui a matrix de ocorrencia pelos valores de div total do Freck
  # mean freck
  
  ## step 3 - prepare output data
  
  
  #' @return A list with the following objects:
  #' \itemize{
  #'     \item{Jetz_per_spp}{Jetz tip-based metric of diversification for each species}
  #'      \item{Jetz_species_site}{Matrix containing
  #'         the values of Jetz tip-based diversification metric for each assemblage} 
  #'     \item{Jetz_harmonic_mean_site}{Harmonic mean of Jetz tip-based
  #'         diversification metric}
  #'      \item{insitu_Jetz_species_sites}{Matrix with model-based Jetz metric for each species
  #'         in each assemblage} 
  #'     \item{insitu_Jetz_harmonic_mean_site}{Harmonic mean of model-based Jetz diversification}
  #' }
  
  if(all(diversification == c("jetz", "freckleton")) == TRUE){
    list_res$Freckleton <- list_res_freckleton
    list_res$Jetz <- list_res
    return(list_res)
  }
  
  if(diversification == "jetz"){
    return(list_res_jetz)
  }
  
  if(diversification == "freckleton"){
    return(list_res_freckleton)
  }
  
   
}