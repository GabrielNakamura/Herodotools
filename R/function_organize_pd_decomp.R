#' Organize output from PD decomposition 
#'
#' @param list1.decomp data frame with PD components
#' @param list2.treetable tibble with edge partition for each community
#'
#' @return
#' @export
#'
#' @examples
organize_pd_decomp <- 
  function(list1.decomp, list2.treetable){
    df_res <- do.call(rbind, list1.decomp)
    list_res_final <- vector(mode = "list")
    list_res_final$PD_decomposition <- df_res
    list_res_final$tree_table_potential <- list2.treetable
    
    # matrix of PD decomposition - dense format
    PD_decomposition_wide <- tidyr::pivot_wider(data = list_res_final$PD_decomposition, names_from = "partition", values_from = "value")
    list_res_final$dense_matrix_PD <- PD_decomposition_wide
    return(list_res_final)
  }