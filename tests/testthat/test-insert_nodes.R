

# create data ----
set.seed(4523)
tree <- ape::rcoal(5)

# Create node area table
node_area <- data.frame(
  area = c("A", "A", "BC", "D"), 
  row.names = paste0("N", 6:9)
)

# Find edge to insert on (e.g., tip 3)
gdata <- ggtree::ggtree(tree)$data
ins_data <- gdata |> dplyr::filter(node %in% c(3, 8, 7)) 
ins_data <- tibble::add_row(ins_data, ins_data[3,])

# Define insertion
# insert 2 nodes in 2 branches
inserts <- tibble::tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.2, 0.8, 0.5, .9),
  node_area = c("C", "AB", "AD", "D")
)

# Run function
result <- insert_nodes(tree, inserts, node_area = node_area)

tree_out <- result$phylo
node_area_out <- result$node_area



test_that("After node addition, tree has the same total length", {
  expect_equal(
    max(ape::node.depth.edgelength(tree)),
    max(ape::node.depth.edgelength(tree_out))
    )
})

test_that("Node labels match rownames of node_area after insertion", {
# Check node labels match
expect_true(!is.null(tree_out$node.label))
expect_true(setequal(tree_out$node.label, rownames(node_area_out)))
expect_equal(tree_out$node.label, rownames(node_area_out))
})

test_that("Original node labels are preserved when no insertions are made", {
  
  original_labels <- paste("nodes", 1:ape::Nnode(tree))
  tree_nl <- tree
  tree_nl$node.label <- original_labels
  
  node_area_nl <- node_area
  rownames(node_area_nl) <- original_labels
  
  tree_no_insert <- insert_nodes(tree_nl, inserts[0, ], node_area_nl)$phylo
  
  expect_equal(tree_no_insert$node.label[1:length(original_labels)], original_labels)
})

test_that("Invalid rownames in node_area throw an error", {
  bad_node_area <- node_area
  rownames(bad_node_area)[1] <- "WRONG"
  expect_error(insert_nodes(tree, inserts, bad_node_area))
})

test_that("Node count increases by number of insertions", {
  n_original_nodes <- tree$Nnode
  n_new_nodes <- nrow(inserts)
  expect_equal(tree_out$Nnode, n_original_nodes + n_new_nodes)
})

test_that("Insertion with invalid event_time throws error", {
  invalid_inserts <- inserts
  invalid_inserts$event_time[1] <- 0  # too small
  expect_error(insert_nodes(tree, invalid_inserts, node_area))
  
  invalid_inserts$event_time[1] <- 1e6  # too large
  expect_error(insert_nodes(tree, invalid_inserts, node_area))
})


test_that("Inserted nodes are labeled with 'new_N' prefix", {
  ana_nodes <- tree_out$node.label[grepl("^new_N", tree_out$node.label)]
  expect_true(length(ana_nodes) == nrow(node_area_out) - tree$Nnode)
})

test_that("Invalid rownames in node_area throw an error", {
  bad_node_area <- node_area
  rownames(bad_node_area)[1] <- "WRONG"
  expect_error(insert_nodes(tree, inserts, bad_node_area))
})


# Branch Length Consistency per edge (inserted vs original)
# This is a bit more elaborate and focuses on per-branch logic


# function to the tests
get_path_between_nodes <- function(gdata, start, end) {
  path <- c(start)
  current <- start
  
  while (current != end) {
    parent <- gdata$parent[gdata$node == current]
    if (is.na(parent)) {
      stop("End node is not an ancestor of the start node.")
    }
    path <- c(path, parent)
    current <- parent
  }
  
  final_path <- path[-length(path)]
  # Return rows in gdata for the nodes along the path
  gdata %>% dplyr::filter(node %in% final_path)
}

test_that("Branch length is preserved per edge after insertion", {
  gdata_orig <- ggtree::ggtree(tree)$data
  gdata_new <- ggtree::ggtree(tree_out)$data
  
  # Unique parent-child pairs (edges) where insertions occurred
  branch_keys <- inserts %>% dplyr::distinct(parent, child)
  
  for (i in seq_len(nrow(branch_keys))) {
    parent_i <- branch_keys$parent[i]
    child_i <- branch_keys$child[i]
    
    # Get original branch length
    bl_orig <- gdata_orig$branch.length[gdata_orig$node == child_i]
    
    # Reconstruct the full chain from child to parent
    path_df <- get_path_between_nodes(gdata_new, start = child_i, end = parent_i)
    
    total_bl <- sum(path_df$branch.length)
    
    expect_equal(round(total_bl, 10), round(bl_orig, 10))
  }
})




test_that("Multiple insertions on same branch are in increasing x order", {
  gdata <- ggtree::ggtree(tree_out)$data
  
  branch_keys <- inserts %>% dplyr::distinct(parent, child)
  
  for (i in seq_len(nrow(branch_keys))) {
    parent_i <- branch_keys$parent[i]
    child_i <- branch_keys$child[i]
    
    # Get full path of nodes from child up to parent
    path_df <- get_path_between_nodes(gdata, start = child_i, end = parent_i)
    
    # Only consider inserted (anagenetic) nodes in the path
    ana_nodes <- path_df %>% 
      dplyr::filter(grepl("^ana_N", label)) %>% 
      dplyr::arrange(x)  # should already be in increasing order of x
    
    # Confirm increasing x
    if (nrow(ana_nodes) > 1) {
      expect_true(all(diff(ana_nodes$x) > 0))
    }
  }
})

test_that("Inserted nodes are at the correct x (time) position", {
  gdata <- ggtree::ggtree(tree_out)$data
  
  inserts_sorted <- inserts %>% dplyr::arrange(child, event_time)
  n_inserts <- nrow(inserts_sorted)
  
  for (i in seq_len(n_inserts)) {
    ins <- inserts_sorted[i, ]
    
    # Match expected node number and label
    expected_node <- ape::Ntip(tree) + tree$Nnode + i
    expected_label <- paste0("new_N", expected_node)
    
    inserted_row <- gdata %>% dplyr::filter(label == expected_label)
    expect_equal(nrow(inserted_row), 1)
    
    parent_x <- gdata$x[gdata$node == ins$parent]
    expected_x <- parent_x + ins$event_time
    actual_x <- inserted_row$x
    
    expect_equal(round(actual_x, 10), round(expected_x, 10))
  }
})

## tests for the internal checks in the function

test_that("Fails when node_area row count does not match number of internal nodes", {
  bad_node_area <- node_area[-1, , drop = FALSE]  # Remove a row
  expect_error(
    insert_nodes(tree, inserts, node_area = bad_node_area),
    "Node count mismatch"
  )
})

test_that("Fails when node labels and node_area rownames do not match", {
  bad_node_area <- node_area
  rownames(bad_node_area)[1] <- "WrongLabel"
  expect_error(
    insert_nodes(tree, inserts, node_area = bad_node_area),
    "Mismatch between tree node labels and rownames in node_area"
  )
})

test_that("Fails if event_time is outside valid range", {
  bad_inserts <- inserts
  bad_inserts$event_time[1] <- 5  # Too long
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "event_time must be > 0 and < original branch length"
  )
})

test_that("Fails if child node in inserts does not exist in tree", {
  bad_inserts <- inserts
  bad_inserts$child[1] <- 999  # Nonexistent node
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "do not exist in the tree"
  )
})

test_that("Fails if parent node in inserts does not exist in tree", {
  bad_inserts <- inserts
  bad_inserts$parent[1] <- 999
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "do not exist in the tree"
  )
})

test_that("Fails if parent-child relation in inserts is not in the tree", {
  bad_inserts <- inserts
  bad_inserts$parent[1] <- bad_inserts$parent[1] + 1  # Break the edge
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "refer to invalid parent-child relationships"
  )
})

test_that("Fails if duplicate insertions on the same edge at the same time", {
  # Create a copy of the inserts data with two insertions at the same time on the same branch
  duplicated_edge <- inserts[1, ]
  duplicated_edge$node_area <- "Dummy"  # make it different in content
  
  bad_inserts <- dplyr::bind_rows(inserts, duplicated_edge)
  
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "Duplicated insertion times for the same branch are not allowed"
  )
})


test_that("Fails if duplicate insertions on the same edge at the same time", {
  bad_inserts <- dplyr::bind_rows(inserts, inserts[1, ])  # Exact duplicate
  expect_error(
    insert_nodes(tree, bad_inserts, node_area = node_area),
    "Duplicated insertion times for the same branch are not allowed"
  )
})

