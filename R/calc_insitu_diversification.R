#' Tip-based in-situ diversification metrics
#'
#' @param W A species occurrence (assemblage) matrix. Rows represent assemblages
#'      (sites) and columns represent species.
#' @param tree A phylogenetic tree object (class `phylo`).
#' @param ancestral.area A one-column data frame indicating the ecorregion of 
#'      occurrence of each node (rows).
#' @param biogeo A data frame with one column indicating the ecoregion of each 
#'      assemblage (row).
#' @param type Character string indicating the type of evolutionary distinctiveness (ED) metric to be used.
#'   Options are `"equal.splits"` (default) or `"fair.proportion"`.
#'
#' @details
#' This function calculates three related diversification metrics:
#' \itemize{
#'   \item **Jetz diversification rate (DR_jetz):** The inverse of the evolutionary distinctiveness (ED) metric
#'         for each species across the full phylogeny. This is the same as 
#'         described in [Jetz et al (2012)](https://www.nature.com/articles/nature11631)
#'         
#'         \deqn{
#'            DR_{i} = (\sum_{j=1}^{N_i}\cdot l_j \cdot 1/2^{j-1})^{-1}
#'         }
#'         
#'   \item **In-situ diversification rate (DR_insitu):** In situ measure of diversification calculated as
#'         the inverse of ED, but restricted to the portion of a species' branches that diversified in the
#'         region (informed in biogeo argument) matching the assemblage.
#'   \item **Proportional diversification rate (DR_prop):** The fraction of diversification that occurred
#'         in a given ecoregion (informed in biogeo argument), calculated as the 
#'         proportion of ED in that region relative to the total ED,
#'         multiplied by the species’ total diversification rate (DR_jetz).
#' }
#'
#' For assemblage-level summaries, site-by-species matrices are generated for each metric, and
#' site-level means are calculated as:
#' \itemize{
#'   \item **Harmonic mean (Jetz and in-situ):** computed across species in each assemblage.
#'     Zeros in the in-situ matrix are ignored in the harmonic mean calculation to
#'     account only for species with diversification in situ.
#'   \item **Arithmetic mean (proportional rates):** used for proportional diversification rates.
#' }
#'
#' @return A named list containing the following elements (if available):
#' \itemize{
#'   \item `jetz_site_sp`: Site-by-species matrix of DR diversification rates.
#'   \item `jetz_comm_mean`: Harmonic mean of DR diversification rates per site.
#'   \item `insitu_site_sp`: Site-by-species matrix of \eqn{DR_in-situ} diversification rates.
#'   \item `insitu_comm_mean`: Harmonic mean of \eqn{DR_in-situ} diversification rates per site
#'         (ignoring zeros when specified).
#'   \item `prop_site_sp`: Site-by-species matrix of proportional diversification rates.
#'   \item `prop_comm_mean`: Arithmetic mean of proportional diversification rates per site.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with a random tree and toy data
#' library(ape)
#' tree <- rcoal(5)
#' W <- matrix(c(1,0,1,1,1,0,1,1,0,1 ), nrow=2, byrow=TRUE,
#'             dimnames = list(c("site1","site2"), tree$tip.label))
#' biogeo <- data.frame(region=c("A","B"))
#' ancestral.area <- data.frame(region=rep("A", tree$Nnode))
#'
#' res <- calc_insitu_diversification(W, tree, ancestral.area, biogeo)
#' str(res)
#' }
#' 
#' @export

calc_insitu_diversification <- function(W,
                                        tree,
                                        ancestral.area, 
                                        biogeo, 
                                        type = "equal.splits"){
  
  if(!is.data.frame(W)) W <- as.matrix(W)
  
  
  if (!is.null(ancestral.area)) {
    if (!is.data.frame(ancestral.area)) {
      stop("ancestral.area must be a data.frame")
    }
    if (ncol(ancestral.area) != 1) {
      stop("ancestral.area must be a data.frame with a single column")
    }
    
    if(!is.null(tree$node.label)){
      if (!all(rownames(ancestral.area) %in% tree$node.label)) {
        stop("Row names of ancestral.area must match node labels in the tree")
      }
    }
    
  }
  
  if(all(type != c("equal.splits", "fair.proportion")) == TRUE){
    stop('"type" must be "equal.splits" or "fair.proportion"')
  }
  
  # step 1 - calc tip rates for all species in all areas ----
  
  # define the unique areas 
  unique_ranges_colapsed <- paste(unique(biogeo[,1]), collapse = "")
  unique_areas <- unique(unlist(strsplit(unique_ranges_colapsed, "")))
  names(unique_areas) <- unique_areas
  
  
  ## calculation for DR_jetz ----
  
  
  ed_total <- calc_ed(tree, type = type)
  jetz_sp <- 1/ed_total
  
  ## calculation for DR_insitu ----
  
  l_ed_insitu <- lapply(unique_areas, function(area){
    calc_ed(tree, ancestral.area, current.area = area, type = type)
  })
  
  l_insitu_sp <- lapply(l_ed_insitu, function(ed_insitu){
    
    in_situ <- 1/ed_insitu
    
    # if ed is zero, DR is zero.
    # meaning the most recent ancestor or range shift
    # is different than the current area of the site. 
    # No diversification in situ!
    in_situ_res <- ifelse(ed_insitu == 0, 0, in_situ)
    
    in_situ_res
  })
  
  ## calculation for DR_prop ----
  
  l_ed_prop <- lapply(l_ed_insitu, function(ed_insitu){
    ed_insitu/ed_total
  })
  
  l_prop_sp <- lapply(l_ed_prop, function(ed_prop) jetz_sp*ed_prop)
  
  # step 2 - populate a comm matrix (site x sp) with tip rates values ----
  
  # generic initial matrix 
  init_matrix <- matrix(NA,
                        nrow = nrow(W),
                        ncol = ncol(W), 
                        dimnames = dimnames(W)
  )
  
  ## populate jetz_site_sp ----
  
  # Initialize an empty matrix with the same dimensions and names as W
  jetz_site_sp <- init_matrix
  
  for(i in colnames(W)){
    # input Jetz total diversification values in occurence matrix
    jetz_site_sp[ , i] <- ifelse(as.numeric(W[ , i]) == 1, jetz_sp[i], NA) 
  }
  
  ## populate insitu_site_sp and prop_site_sp ----
  insitu_site_sp <- init_matrix
  prop_site_sp <- init_matrix
  
  # Loop through each community
  for (i in seq_len(nrow(W))) {
    region <- biogeo[i, 1]
    
    # Populate insitu_site_sp if needed 
    if (exists("l_insitu_sp")) {
      insitu_vals <- l_insitu_sp[[region]][colnames(W)]
      insitu_site_sp[i, ] <- ifelse(as.numeric(W[i, ]) == 1, as.numeric(W[i, ]) * insitu_vals, NA)
      
    }
    
    # Populate prop_site_sp if needed
    if (exists("l_prop_sp")) {
      prop_vals <- l_prop_sp[[region]][colnames(W)]
      prop_site_sp[i, ] <- ifelse(
        as.numeric(W[i, ]) == 1, as.numeric(W[i, ]) * prop_vals, NA
      )
      
      
    }
  }
  
  # Step 3 - compute community mean  -----
  
  ## jetz_comm_mean ----
  
  jetz_comm_mean <- apply(jetz_site_sp, 1,  function(x){
    calc_harmonic_mean(x, na.rm = T)
  })
  
  ## insitu_comm_mean ----
  insitu_comm_mean <- apply(insitu_site_sp, 1, function(x){
    
    # By ignoring zero the site harmonic mean takes into account only the 
    # species with in situ diversification rates. 
    calc_harmonic_mean(x, na.rm = T, ignore.zero = T)
  })
  
  ## prop_comm_mean ----
  sum_prop_site <- rowSums(prop_site_sp, na.rm = T)
  sum_jetz_site <- rowSums(jetz_site_sp, na.rm = T)
  prop_comm_mean <- sum_prop_site/sum_jetz_site
  
  
  # Step 4 - prepare results list
  
  list_res <- list()
  
  # Add elements only if they exist
  if (exists("jetz_site_sp"))   list_res$jetz_site_sp   <- jetz_site_sp
  if (exists("jetz_comm_mean")) list_res$jetz_comm_mean <- jetz_comm_mean
  
  if (exists("insitu_site_sp"))   list_res$insitu_site_sp   <- insitu_site_sp
  if (exists("insitu_comm_mean")) list_res$insitu_comm_mean <- insitu_comm_mean
  
  if (exists("prop_site_sp"))   list_res$prop_site_sp   <- prop_site_sp
  if (exists("prop_comm_mean")) list_res$prop_comm_mean <- prop_comm_mean
  
  return(list_res)
  
}