#' Calculate a evolutionary distinctiveness (ED) based on biogeographical model
#' 
#' For each tip, the function calculates the ED based on the current area and the 
#' ancestral area based on a biogeographical model. ED is calculated for the connecting 
#' branched from tip to internal node that are in the same area as the current area. 
#' This means that the evolutionary distinctiveness occurred in situ (same biogeographical region). 
#' If current area is not provided, the function gives same result as 'picante::evol.distinct()'
#'
#' @param tree Phylogenetic tree of class 'phylo'.
#' @param ancestral.area One column data frame indicating the Ecoregion of occurrence of each node (rows)
#' @param current.area character. A single area that will represent the current area of all tips in the phylogeny.
#' @param type Character indicating the type of evolutionary distinctiveness metric to be used. 
#'    It can be "equal.splits" or "fair.proportion", "equal-splits" is the default.
#'
#' @returns a named numeric vector with values for all tips in the phylogeny
#' @export
#'
#' @examples
#' 
#' # example for calc_ed
#' # generate simlutated data ----
#' set.seed(4523)
#' tree_sim <- ape::rcoal(5)
#'
#' # Create node area table (only bifurcating nodes)
#' node_area <- data.frame(
#'     area = c("A", "A", "BC", "AD"),
#'     row.names = paste0("N", 6:9)
#'     )
#'     
#' # create insertion of nodes to simulate for change in area in the branch 
#' # Use ggtree to identify insertion points 
#' 
#' gdata <- ggtree::ggtree(tree_sim)$data
#' ins_data <- gdata %>% dplyr::filter(node %in% c(3, 9, 8))
#' ins_data <- tibble::add_row(ins_data, ins_data[2, ]) %>% arrange(node)
#' 
#' # Define inserted nodes with areas
#' inserts <- tibble::tibble(
#'   parent = ins_data$parent,
#'   child = ins_data$node,
#'   event_time = c(0.2, 0.5, 0.9, .3),
#'   node_area = c("AC", "AB", "BC", "D")
#'   )
#'   
#' # Insert nodes 
#' result <- insert_nodes(tree_sim, inserts, node_area = node_area)
#' tree_out <- result$phylo
#' node_area_out <- result$node_area
#' 
#' 
#' # Compute Evolutionary Distinctiveness --------
#' 
#' # ED total (same results as picante::evol.distinct())
#' ed_total <- calc_ed(  
#'   tree_out, 
#'   node_area = node_area_out, 
#'   type = "equal.splits")
#'   
#' # ED partial. Account for in situ distinctiveness only
#' 
#' ed_partial_A <- calc_ed(
#'   tree_out,
#'   node_area = node_area_out, 
#'   current_area = "A", 
#'   type = "equal.splits")
#'   
#' ed_partial_D <- calc_ed(
#'   tree_out,
#'   node_area = node_area_out, 
#'   current_area = "D", 
#'   type = "equal.splits")
#'   
#' data.frame(
#'   species = names(ed_total), 
#'   ed_total,
#'   ed_partial_A,
#'   ed_partial_D
#'   )
#'   
#' # tree visualization
#' tree_viz <- tree_out
#' tree_viz$node.label <- node_area_out$area
#' 
#' plot(tree_viz, show.node.label = T)
#' nodelabels(adj = 1.2)
#' 
#' 
#' 


calc_ed <- function(tree, ancestral.area, current.area = NULL, type = c("equal.splits", "fair.proportion")) {
  type <- match.arg(type)
  
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  if (!is.data.frame(ancestral.area)) stop("ancestral.area must be a data.frame with rownames as node labels and a column 'area'")
  
  edge <- tree$edge
  edge.length <- tree$edge.length
  n_tips <- length(tree$tip.label)
  n_node <- ape::Nnode(tree)
  
  # Ensure node labels are present on the tree
  if (is.null(tree$node.label)) {
    tree$node.label <- paste0("N", seq_len(tree$Nnode) + ape::Ntip(tree))
  }
  
  # Internal node label check ---
  expected_nums <- (length(tree$tip.label) + 1):(length(tree$tip.label) + ape::Nnode(tree))
  expected_labels <- c(
    paste0("N", expected_nums),
    paste0("new_N", expected_nums)
  )
  
  invalid_labels <- which(!(tree$node.label %in% expected_labels))
  if (length(invalid_labels) > 0) {
    stop("Invalid node labels found: ", paste(tree$node.label[invalid_labels], collapse = ", "))
  }
  
  
  # Map node labels to node numbers
  node_labels <- setNames(n_tips + seq_along(tree$node.label), tree$node.label)
  
  # Determine ignored nodes
  ignore_nodes <- grep("^new_N", tree$node.label, value = TRUE)
  ignored_node_nums <- unname(node_labels[ignore_nodes])
  
  # Convert ancestral.area to map from node number to area
  ancestral.area$node <- node_labels[rownames(ancestral.area)]
  area_by_node <- setNames(as.character(ancestral.area$area), ancestral.area$node)
  
  # For fair proportion: precompute descendant tip counts
  if (type == "fair.proportion") {
    descendant_tips <- integer(n_tips + n_node)
    for (node in 1:(n_tips + n_node)) {
      descendant_tips[node] <- length(phangorn::Descendants(tree, node, "tips")[[1]])
    }
  }
  
  ed <- setNames(numeric(n_tips), tree$tip.label)
  
  for (tip in 1:n_tips) {
    node <- tip
    ed_value <- 0
    decay_index <- 0
    pending_length <- 0
    
    # walk at nodes from tip to root while sum the ed values
    while (node != (n_tips + 1)) {
      parent_edge_index <- which(edge[, 2] == node)
      parent_node <- edge[parent_edge_index, 1]
      
      # check if node area is different than the site area (current.area)
      if (!is.null(current.area)) {
        parent_area <- area_by_node[as.character(parent_node)]
        if (is.na(parent_area) || !grepl(paste0("[", current.area, "]"), parent_area)) {
          if (type == "equal.splits") {
            ed_value <- ed_value + (pending_length * (0.5 ^ decay_index))
          }
          if (type == "fair.proportion") {
            ed_value <- ed_value + (pending_length / descendant_tips[node])
          }
          
          break
        }
      }
      
      branch_len <- edge.length[parent_edge_index]
      pending_length <- pending_length + branch_len
      
      # calculate the ed value 
      if (!(parent_node %in% ignored_node_nums)) {
        if (type == "equal.splits") {
          ed_value <- ed_value + (pending_length * (0.5 ^ decay_index))
          pending_length <- 0
          decay_index <- decay_index + 1
        } else if (type == "fair.proportion") {
          if (node <= n_tips) {
            ed_value <- ed_value + pending_length
          } else {
            ed_value <- ed_value + (pending_length / descendant_tips[node])
          }
          pending_length <- 0
        }
      }
      
      # go next node
      node <- parent_node
    }
    
    ed[tip] <- ed_value
  }
  
  return(ed)
}
