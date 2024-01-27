
#' Internal function used to calculate PD components
#'
#' @param data a tibble containing the phylogenetic object 
#' @param comm community matrix
#' @param phy newick phylogenetic object
#'
#' @return
#' @export
#'
#' @examples
#' 
partition_pd <- 
  function(data, comm, phy){
    comm_names <- rownames(comm)
    
    data2 <- 
      data %>% 
      dplyr::mutate(pres = ifelse(label %in% colnames(comm), "pres", "abs")) %>% 
      dplyr::mutate(ancestor = data$label[data$parent]) %>% 
      dplyr::mutate(descendant = data$label[data$node]) %>% 
      dplyr::mutate(ancestor1 = gsub(pattern = ".*_", replacement = "", x = ancestor)) %>% 
      dplyr::mutate(descendent1 = gsub(pattern = ".*_", replacement = "", x = descendant)) %>% 
      dplyr::mutate(partition.IS = ifelse(ancestor1 == "IS" & descendent1 == "IS" |
                                            ancestor1 == "IS" & descendent1 == "EM" |
                                            ancestor1 == "EM" & descendent1 == "IS" | 
                                            ancestor1 == "EM" & descendent1 == "EM" |
                                            ancestor1 == "IS" & pres == "pres" |
                                            ancestor1 == "EM" & pres == "pres",
                                          "IS", NA)) %>% 
      dplyr::mutate(partition.IM = ifelse(ancestor1 == "IM" & descendent1 == "IS" |
                                            ancestor1 == "IM" & descendent1 == "EM" |
                                            ancestor1 == "ESD" & descendent1 == "IS" |
                                            ancestor1 == "ESD" & descendent1 == "EM"|
                                            ancestor1 == "ESD" & pres == "pres"| 
                                            ancestor1 == "IM" & pres == "pres",
                                          "IM", NA)) %>%
      dplyr::mutate(partition.EM = ifelse(ancestor1 == "IS" & descendent1 == "IM" |
                                            ancestor1 == "EM" & descendent1 == "ESD" |
                                            ancestor1 == "EM" & descendent1 == "IM" | 
                                            ancestor1 == "IS" & descendent1 == "ESD" |
                                            ancestor1 == "IS" & pres == "abs" & descendent1 != "IS" & descendent1 != "IM" & descendent1 != "ESD" & descendent1 != "EM" | 
                                            ancestor1 == "EM" & pres == "abs" & descendent1 != "IS" & descendent1 != "IM" & descendent1 != "ESD" & descendent1 != "EM", 
                                          "EM", NA)) %>% 
      dplyr::mutate(partition.ESD = ifelse(ancestor1 == "ESD" & descendent1 == "ESD" |
                                             ancestor1 == "ESD" & descendent1 == "IM" |
                                             ancestor1 == "IM" & descendent1 == "IM" |
                                             ancestor1 == "IM" & descendent1 == "ESD" |
                                             ancestor1 == "IM" & pres == "abs" & descendent1 != "IS" & descendent1 != "IM" & descendent1 != "ESD" & descendent1 != "EM"| 
                                             ancestor1 == "ESD" & pres == "abs" & descendent1 != "IS" & descendent1 != "IM" & descendent1 != "ESD" & descendent1 != "EM", 
                                           "ESD", NA))
    
    # adding a group to all species - this will be useful to plot the partitions
    data2 <- 
      data2 %>% 
      mutate(group = coalesce(partition.IS, partition.IM, partition.EM, partition.ESD))  
    
    # calculating PD components
    PDinsitu <- 
      data2 %>% 
      dplyr::filter(partition.IS == "IS") %>% 
      dplyr::select(branch.length) %>% 
      sum(na.rm = T)
    
    PDimmigration <- 
      data2 %>% 
      dplyr::filter(partition.IM == "IM") %>% 
      dplyr::select(branch.length) %>% 
      sum(na.rm = T)
    
    PDemigration <- 
      data2 %>% 
      dplyr::filter(partition.EM == "EM") %>% 
      dplyr::select(branch.length) %>% 
      sum(na.rm = TRUE)
    
    PDexsitu <- 
      data2 %>% 
      dplyr::filter(partition.ESD == "ESD") %>% 
      dplyr::select(branch.length) %>% 
      sum(na.rm = TRUE)
    
    
    PDtotal <- PDinsitu + PDimmigration + PDemigration + PDexsitu # adding up all components
    PDfaith <- picante::pd(samp = comm, tree = phy, include.root = FALSE)[1, 1]  # including PD calculated accordingly to Faith approach
    
    #joining all results
    data_res <- 
      data.frame(partition = c("PDinsitu", "PDimmigration", "PDemigration", "PDexsitu", "PDtotal", "PDfaith"), 
                 value = c(PDinsitu, PDimmigration, PDemigration, PDexsitu, PDtotal, PDfaith), community = comm_names)
    
    
    res_list <- vector(mode = "list")
    res_list$data_res <- data_res
    res_list$table_tree <- data2
    return(res_list)
  }