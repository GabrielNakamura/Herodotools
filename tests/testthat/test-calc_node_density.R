library(testthat)
library(ape)
library(dplyr)
library(tibble)
library(ggtree)
library(picante)

# Set up shared test data
set.seed(4523)
tree_sim <- ape::rcoal(5)

# Create node area table (only bifurcating nodes)
node_area <- data.frame(
  area = c("A", "A", "BC", "D"),
  row.names = paste0("N", 6:9)
)

# Use ggtree to identify insertion points
gdata <- ggtree::ggtree(tree_sim)$data
ins_data <- gdata %>% dplyr::filter(node %in% c(3, 9, 8))
ins_data <- tibble::add_row(ins_data, ins_data[2, ]) %>% arrange(node)

# Define inserted nodes with areas
inserts <- tibble::tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.2, 0.5, 0.9, .3),
  node_area = c("AC", "AB", "BC", "AD")
)

# Insert nodes
result <- insert_nodes(tree_sim, inserts, node_area = node_area)
tree_out <- result$phylo
node_area_out <- result$node_area

plot(tree_out)
nodelabels()


node_density_metric <- function(tree) {
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  
  edge <- tree$edge
  n_tips <- length(tree$tip.label)
  
  # Compute node depths from root
  node_depths <- ape::node.depth.edgelength(tree)
  
  # Initialize output
  nd <- setNames(numeric(n_tips), tree$tip.label)
  
  for (tip in 1:n_tips) {
    node <- tip
    node_count <- 0
    
    while (node != (n_tips + 1)) {
      parent_edge_index <- which(edge[, 2] == node)
      parent_node <- edge[parent_edge_index, 1]
      
      node_count <- node_count + 1
      node <- parent_node
    }
    
    # Use depth of tip as denominator (age from root)
    nd[tip] <- node_count / node_depths[tip]
  }
  

  return(nd)
}


# tests


test_that("Node density calculation is correct when using the original formulation", {
  nd_custom <- calc_node_density(tree_sim, ancestral.area = node_area)
  nd_original <- node_density_metric(tree_sim)
  expect_equal(nd_custom, nd_original)
})


test_that("Node density computes contributions before stopping", {
  el_df <- data.frame(
    tree_out$edge,
    edge.length = tree_out$edge.length,
    labels = c(tree_out$tip.label, tree_out$node.label[-1])
  )
  
  ed_partialB <- calc_node_density(tree_out,
                                   ancestral.area = node_area_out, 
                                   current.area = "B")
  ed_partialC <- calc_node_density(tree_out,
                                   ancestral.area = node_area_out,
                                   current.area = "C")
    
    edge_1 <- el_df %>%
      filter(labels == "t3") %>%
      pull(edge.length)
    edge_2 <- el_df %>%
      filter(labels == "N8") %>%
      pull(edge.length)
    edge_3 <- el_df %>%
      filter(labels == "new_N12") %>%
      pull(edge.length)
    
    
    # nd = n_nodes / sum of edge_length
    ed_t3_B <- 1/(edge_1 + edge_2 + edge_3)
    ed_t3_C <- 1/(edge_1 + edge_2)
    
    expect_equal(as.numeric(ed_partialB["t3"]), ed_t3_B)
    expect_equal(as.numeric(ed_partialC["t3"]), ed_t3_C)
  
})

test_that("Node density is zero when current.area is differnt from the 
ancestral area of the most recent ancestor", {
  
  ed_partialC <- calc_node_density(tree_out,
                                   ancestral.area = node_area_out,
                                   current.area = "C")
  
  
  expect_equal(as.numeric(ed_partialC["t2"]), 0)        
  expect_equal(as.numeric(ed_partialC["t4"]), 0)          

  })


test_that("Node density is zero when last area shit occured after the most
recent speciation event", {
  
  ed_partialC <- calc_node_density(tree_out,
                                   ancestral.area = node_area_out,
                                   current.area = "C")
  
  
  expect_equal(as.numeric(ed_partialC["t1"]), 0)          
  
})





