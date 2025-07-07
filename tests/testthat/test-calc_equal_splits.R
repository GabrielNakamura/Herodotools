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
ins_data <- tibble::add_row(ins_data, ins_data[2,]) %>% arrange(node)

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

test_that("Total ED matches picante::evol.distinct()", {
  ed_custom <- calc_equal_splits_ed(tree_out, node_area = node_area_out)
  ed_picante <- picante::evol.distinct(tree_sim, type = "equal.splits")
  expect_equal(ed_custom, setNames(ed_picante$w, ed_picante$Species))
})

test_that("Partial ED is less than total ED when path is interrupted", {
  ed_total <- calc_equal_splits_ed(tree_out, node_area = node_area_out)
  ed_partial <- calc_equal_splits_ed(tree_out, node_area = node_area_out, current_area = "B")
  
  # Tip 3 should stop at node with area "AD"
  expect_true(ed_partial["t3"] < ed_total["t3"])
  
  # Tip 1 should stop at node with area "AB"
  expect_true(ed_partial["t1"] < ed_total["t1"])
})

test_that("Partial ED equals total ED when no interruption occurs", {
  ed_total <- calc_equal_splits_ed(tree_out, node_area = node_area_out)
  ed_partial <- calc_equal_splits_ed(tree_out, node_area = node_area_out, current_area = "A")
  
  # Assuming tip 4 has no inserted nodes with area incompatible with "B"
  expect_equal(ed_partial["t1"], ed_total["t1"])
})

test_that("Partial ED computes contributions before stopping", {
  el_df <- data.frame(
    tree_out$edge,
    edge.length = tree_out$edge.length, 
    labels = c(tree_out$tip.label, tree_out$node.label[-1])
  )
  
  ed_partialB <- calc_equal_splits_ed(tree_out, node_area = node_area_out, current_area = "B")
  ed_partialC <- calc_equal_splits_ed(tree_out, node_area = node_area_out, current_area = "C")
  
  # Compute expected ED for t3, manually but without hardcoding indices
  # Step 1: edge leading to t3
  edge_1 <- el_df %>% 
    filter(labels == "t3") %>% 
    pull(edge.length)
  
  # Step 2: edge from t3's parent to bifurcating node with area containing B or C
  edge_2 <- el_df %>% 
    filter(labels == "N8") %>% 
    pull(edge.length)
  
  # Step 3: additional edge if area == B
  edge_3 <- el_df %>%
    filter(labels == "new_N12") %>%
    pull(edge.length)
  
  # Expected ED values (equal splits weights applied)
  ed_t3_B <- edge_1 + ((edge_2 + edge_3) * 0.5)
  ed_t3_C <- edge_1 + (edge_2 * 0.5)
  
  # Ensure partial ED is non-zero for tips that get interrupted mid-path
  expect_equal(as.numeric(ed_partialB["t3"]), ed_t3_B)
  expect_equal(as.numeric(ed_partialC["t3"]), ed_t3_C)
})


