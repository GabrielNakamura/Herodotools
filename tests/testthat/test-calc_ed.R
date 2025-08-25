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


# tests
test_that("calc_ed requires both or none of ancestral.area and current.area", {
  # Case 1: Only ancestral.area provided → should fail
  expect_error(
    calc_ed(tree_sim, ancestral.area = node_area, type = "equal.splits"),
    "both 'ancestral.area' and 'current.area'"
  )
  
  # Case 2: Only current.area provided → should fail
  expect_error(
    calc_ed(tree_sim, current.area = "A", type = "equal.splits"),
    "both 'ancestral.area' and 'current.area'"
  )
  
  # Case 3: Both provided → should work (no error)
  expect_silent(
    calc_ed(tree_sim, ancestral.area = node_area, current.area = "A", type = "equal.splits")
  )
  
  # Case 4: None provided → should work (no error)
  expect_silent(
    calc_ed(tree_sim, type = "equal.splits")
  )
})



test_that("Total ED matches picante::evol.distinct(), without insert nodes ", {
  ed_custom <- calc_ed(tree_sim, type = "equal.splits")
  ed_picante <- picante::evol.distinct(tree_sim, type = "equal.splits")
  expect_equal(ed_custom, setNames(ed_picante$w, ed_picante$Species))
  
  
  ed_custom <- calc_ed(tree_sim, type = "fair.proportion")
  ed_picante <- picante::evol.distinct(tree_sim, type = "fair.proportion")
  expect_equal(ed_custom, setNames(ed_picante$w, ed_picante$Species))
})


test_that("Total ED matches picante::evol.distinct(), with inserted nodes ", {
  ed_custom <- calc_ed(tree_out, type = "equal.splits")
  ed_picante <- picante::evol.distinct(tree_sim, type = "equal.splits")
  expect_equal(ed_custom, setNames(ed_picante$w, ed_picante$Species))
  
  ed_custom <- calc_ed(tree_out, type = "fair.proportion")
  ed_picante <- picante::evol.distinct(tree_sim, type = "fair.proportion")
  expect_equal(ed_custom, setNames(ed_picante$w, ed_picante$Species))
})


test_that("Partial ED is less than total ED when path is interrupted", {
  for (type in c("equal.splits", "fair.proportion")) {
    ed_total <- calc_ed(tree_out, type = type)
    ed_partial <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "B", type = type)
    
    expect_true(ed_partial["t3"] < ed_total["t3"], info = paste("Type:", type, "Tip: t3"))
    expect_true(ed_partial["t1"] < ed_total["t1"], info = paste("Type:", type, "Tip: t1"))
  }
})

test_that("Partial ED equals total ED when no interruption occurs", {
  for (type in c("equal.splits", "fair.proportion")) {
    ed_total <- calc_ed(tree_out, type = type)
    ed_partial <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "A", type = type)
    
    expect_equal(ed_partial["t1"], ed_total["t1"], info = paste("Type:", type, "Tip: t1"))
  }
})

test_that("Partial ED (equal.splits) computes contributions before stopping", {
  el_df <- data.frame(
    tree_out$edge,
    edge.length = tree_out$edge.length,
    labels = c(tree_out$tip.label, tree_out$node.label[-1])
  )
  
  for (type in c("equal.splits", "fair.proportion")) {
    ed_partialB <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "B", type = "equal.splits")
    ed_partialC <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "C", type = "equal.splits")
    
    edge_1 <- el_df %>%
      filter(labels == "t3") %>%
      pull(edge.length)
    edge_2 <- el_df %>%
      filter(labels == "N8") %>%
      pull(edge.length)
    edge_3 <- el_df %>%
      filter(labels == "new_N12") %>%
      pull(edge.length)
    
    ed_t3_B <- edge_1 + ((edge_2 + edge_3) * 0.5)
    ed_t3_C <- edge_1 + (edge_2 * 0.5)
    
    expect_equal(as.numeric(ed_partialB["t3"]), ed_t3_B)
    expect_equal(as.numeric(ed_partialC["t3"]), ed_t3_C)
  }
})


test_that("Partial ED (fair.proportion) computes contributions before stopping", {
  el_df <- data.frame(
    tree_out$edge,
    edge.length = tree_out$edge.length,
    labels = c(tree_out$tip.label, tree_out$node.label[-1])
  )
  
  ed_partialB <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "B", type = "fair.proportion")
  ed_partialC <- calc_ed(tree_out, ancestral.area = node_area_out, current.area = "C", type = "fair.proportion")
  
  # Step 1: edge leading to t3 (tip edge, full length)
  edge_1 <- el_df %>%
    filter(labels == "t3") %>%
    pull(edge.length)
  
  # Step 2: edge from t3's parent to N8 (internal edge)
  edge_2 <- el_df %>%
    filter(labels == "N8") %>%
    pull(edge.length)
  desc_N8 <- length(phangorn::Descendants(tree_out, which(tree_out$node.label == "N8") + length(tree_out$tip.label), "tips")[[1]])
  
  # Step 3: edge from N8 to new_N12 (internal edge)
  edge_3 <- el_df %>%
    filter(labels == "new_N12") %>%
    pull(edge.length)
  desc_new_N12 <- length(phangorn::Descendants(tree_out, which(tree_out$node.label == "new_N12") + length(tree_out$tip.label), "tips")[[1]])
  
  # Expected ED values (fair proportion weights applied)
  ed_t3_B <- edge_1 + ((edge_2 + edge_3) / desc_N8)
  ed_t3_C <- edge_1 + (edge_2 / desc_N8) # stops before new_N12
  
  expect_equal(as.numeric(ed_partialB["t3"]), ed_t3_B)
  expect_equal(as.numeric(ed_partialC["t3"]), ed_t3_C)
})


test_that("calc_ed stops with invalid node labels", {
  bad_tree <- tree_sim
  bad_tree$node.label <- c("N6", "N7_extra", "new_N8", "new_N9")
  expect_error(
    calc_ed(bad_tree, type = "equal.splits"),
    regexp = "Invalid node labels"
  )
})


test_that("calc_ed fails if ancestral.area is not a data.frame", {
  tree <- ape::rcoal(5)
  tree$node.label <- paste0("N", seq_len(tree$Nnode) + ape::Ntip(tree))
  
  expect_error(
    calc_ed(tree, ancestral.area = list(a = 1), current.area = "A"),
    "ancestral.area must be a data.frame"
  )
})

test_that("calc_ed fails if ancestral.area has more than one column", {
  tree <- ape::rcoal(5)
  tree$node.label <- paste0("N", seq_len(tree$Nnode) + ape::Ntip(tree))
  
  aa <- data.frame(a = rep("A", tree$Nnode), b = 1:tree$Nnode)
  rownames(aa) <- tree$node.label
  
  expect_error(
    calc_ed(tree, ancestral.area = aa, current.area = "A"),
    "ancestral.area must be a data.frame with a single column"
  )
})

test_that("calc_ed fails if rownames of ancestral.area do not match tree$node.label", {
  tree <- ape::rcoal(5)
  tree$node.label <- paste0("N", seq_len(tree$Nnode) + ape::Ntip(tree))
  
  aa <- data.frame(a = rep("A", tree$Nnode))
  rownames(aa) <- paste0("X", seq_len(tree$Nnode))  # wrong labels
  
  expect_error(
    calc_ed(tree, ancestral.area = aa, current.area = "A"),
    "Row names of ancestral.area must match node labels in the tree"
  )
})

test_that("calc_ed works silently when ancestral.area is valid", {
  tree <- ape::rcoal(5)
  tree$node.label <- paste0("N", seq_len(tree$Nnode) + ape::Ntip(tree))
  
  aa <- data.frame(foo = rep("A", tree$Nnode))
  rownames(aa) <- tree$node.label
  
  expect_silent(
    calc_ed(tree, ancestral.area = aa, current.area = "A")
  )
})

test_that("calc_ed generates node labels when tree$node.label is NULL", {
  tree <- ape::rcoal(5)  # no node labels
  aa <- data.frame(bar = rep("A", tree$Nnode))
  
  expect_silent(
    res <- calc_ed(tree, ancestral.area = aa, current.area = "A")
  )
  
  # After the function, labels should have been auto-generated
  expect_true(all(grepl("^N[0-9]+", tree$node.label)))
})

test_that("calc_ed overwrites ancestral.area rownames when tree$node.label is NULL", {
  tree <- ape::rcoal(5)  # no node labels
  aa <- data.frame(bar = rep("A", tree$Nnode))
  rownames(aa) <- paste0("wrong", seq_len(tree$Nnode))  # mismatched rownames
  
  expect_silent(
    calc_ed(tree, ancestral.area = aa, current.area = "A")
  )
})
