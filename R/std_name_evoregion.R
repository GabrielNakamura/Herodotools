#' Standardize the evoregion names after 'adegenet::find.clusters'
#' 
#' Although the results of 'find.cluster' are consistent when function is runned multiple times, the order of the names of the groups are random. 
#' Here, the function reorder the name of the groups following the sequence the appear in the vector. 
#' 
#'
#' @param clust_evoreg group vector from the 'adegenet::find.clusters'
#' 
#' @return factor vector with the groups
#' 
#' @details
#' The function was designed to be used within the function \code{`calc_evoregions`}
#' 
#' 
#' @examples
#' 
#' @import dplyr
#' @importFrom stringr word 
#' @importFrom stringr fixed 
  
std_name_evoregion <- function(clust_evoreg){
  
  #split groups
  n_clust <- length(levels(clust_evoreg))
  splited_grp <- split(clust_evoreg, clust_evoreg) 
  
  # find the new order
  grp_ord <- lapply(splited_grp, \(x){names(x)[1]}) %>% unlist()
  
  df_grp <- data.frame(
    grp = as.integer(names(grp_ord)), 
    row = as.integer(grp_ord)
  ) %>%
    dplyr::arrange(row) %>% 
    dplyr::mutate(new_grp = 1:n_clust)
  
  # rename the groups
  for(i in 1:nrow(df_grp)){
    grp <- df_grp[i, "grp"]
    new_grp <- df_grp[i, "new_grp"]
    
    splited_grp[[grp]][grp==grp] <- new_grp
  }
  
  # rearrange the vector the sites in the same order of the input data
  unlist_grp <- splited_grp %>% unlist()
  
  nm_vec <- stringr::word(names(unlist_grp), start = 2, sep = stringr::fixed("."))
  
  out <- unlist_grp[order(as.integer(nm_vec))]
  names(out) <- names(clust_evoreg)
  
  # return the result
  out
}