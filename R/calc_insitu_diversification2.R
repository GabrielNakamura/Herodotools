calc_insitu_diversification2 <- function(W,
                                         tree,
                                         ancestral.area, 
                                         biogeo, 
                                         diversification = "jetz",
                                         type = "equal.splits"){
  
  
  if(all(diversification != c("jetz", "freckleton")) == TRUE){
    stop("Tip-based diversification measures must be one of 'jetz' or 'freckleton'")
  }
  
  # # Function to compute harmonic mean including zeros
  # # x is a vector of a site
  # # NA are removed from the calculation 
  # harmonic_mean <- function(x) {
  #   x <- x[!x==0]
  #  
  #   length(x) / sum(1 / x)
  # }
  
 
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
      
      # adds a very low value when the most recent ancestor or range shift
      # is different than the current area of the site
      matrix_XJetz[i, ] <- ifelse(
        W[i, ] == 1 & jetz_vals == 0, 0.00001, W[i, ] * jetz_vals
      )
    }
    
    # Summarazing results for each site ---------------------------------------
    
    ## calculating harmonic mean for Jetz
    
    # sum of communities local diversification
    sum_localDiv <- apply(matrix_XJetz,
                          1,
                          function(x) sum(x[which(x != 0 & x != 0.00001)]))
    
    
    #object to receive Jetz diversification values
    matrix_totalDiv_Jetz <- W
    
    # Total jetz values for each community
    for(i in colnames(W)){
      matrix_totalDiv_Jetz[ , i] <- ifelse(matrix_totalDiv_Jetz[ , i] == 1, 
                                           jetz_total[i],
                                           0) # input Jetz total diversification values in occurence matrix
    }
    
    #sum of communities total diversification
    sum_totalDiv <- apply(matrix_totalDiv_Jetz,
                          1,
                          function(x) sum(x[which(x != 0)]))
    
    numerator_values <- apply(matrix_totalDiv_Jetz, 1, function(x){
      sum(1/x[which(x != 0)])
    })
    
    denom_values <- apply(matrix_totalDiv_Jetz, 1, function(x){
      length(which(x != 0))
    })
    
    #harmonic mean for Jetz total diversification
    JetzTotalComm_harmonic <- denom_values/numerator_values
    
    #harmonic mean for Jetz local diversification
    JetzLocalComm_harmonic <- (JetzTotalComm_harmonic*(sum_localDiv/sum_totalDiv))
    list_res_jetz <- vector(mode = "list", length = 5)
    list_res_jetz[[1]] <- jetz_total
    list_res_jetz[[2]] <- matrix_totalDiv_Jetz
    list_res_jetz[[3]] <- JetzTotalComm_harmonic # jetz harmonic mean diversification per community
    list_res_jetz[[4]] <- matrix_XJetz # model-based metric per community per species
    list_res_jetz[[5]] <- JetzLocalComm_harmonic # model-based metric harmonic mean 
    names(list_res_jetz) <- c("Jetz_per_spp",
                              "Jetz_species_site",
                              "Jetz_harmonic_mean_site",
                              "insitu_Jetz_species_sites", 
                              "insitu_Jetz_harmonic_mean_site")
    }
  

  
  # if freckleton:
  #  # compute total diversification 
  # compute the nd for all the areas as current area 
  
  if(any(diversification == "freckleton")){
  
    nd_total <- calc_node_density(tree, ancestral.area)
    
    l_nd_local <- lapply(unique_areas, function(area){
      calc_node_density(tree, ancestral.area, current.area = area)
      })
    
    matrix_XFreck<- matrix(0,
                           nrow = nrow(W),
                           ncol = ncol(W),
                           dimnames = list(rownames(W), colnames(W)))
    
    # Loop through each community
    for (i in seq_len(nrow(W))) {
      region <- biogeo[i, 1]
      nd_vals <- l_nd_local[[region]][colnames(W)]
      
      # adds a very low value when the most recent ancestor or range shift
      # is different than the current area of the site
      matrix_XFreck[i, ] <- ifelse(
        W[i, ] == 1 & nd_vals == 0, 0.00001, W[i, ] * nd_vals
      )
    }
    
    FreckLocalComm_mean <- sapply(1:nrow(matrix_XFreck), function(i){
      pres <- matrix_XFreck[i, ][W[i, ] >= 1]
      mean(pres, na.rm = T)
    })
    
    #objet to receive Jetz diversification values
    matrix_totalDiv_Freck <- W
    
    #substitui a matrix de ocorrencia pelos valores de div total do Freck
    for(i in colnames(W)){
      matrix_totalDiv_Freck[, i] <- ifelse(matrix_totalDiv_Freck[, i] >= 1, nd_total[i], 0)
    }
    
    FreckTotalComm_mean <- sapply(1:nrow(matrix_totalDiv_Freck), function(i){
      pres <- matrix_totalDiv_Freck[i, ][W[i, ] >= 1]
      mean(pres, na.rm = T)
    })
    
    list_res_freckleton <-  data.frame(FreckTotal = FreckTotalComm_mean,
                                       FreckLocal = FreckLocalComm_mean,
                                       FreckLocalProp = FreckLocalComm_mean/
                                         FreckTotalComm_mean,
                                       row.names = row.names(W))
    
    
    }
  
  
  
 
  
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