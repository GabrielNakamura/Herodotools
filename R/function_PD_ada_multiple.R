#' Calculation of PD using ada reconstruction 
#'
#' @param ada.obj 
#' @param thresold 
#' @param comm 
#'
#' @return
#' @export
#'
#' @examples
PD_ada <- 
  function(ada.obj, thresold, comm){
    comm.name <- rownames(comm)
    res_list <- lapply(comm.name, function(x) calc_PD_ada(ada.obj = ada.obj, threshold = threshold, comm.name = x))
    res_final <- do.call(rbind, res_list)
    return(res_final)
  }

