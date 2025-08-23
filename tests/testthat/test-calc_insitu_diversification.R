# data ----

# Load toy tree ----
data(toy_treeEx)
tree_sim <- toy_treeEx
tree_sim$node.label <- NULL
tree_sim$edge.length <- tree_sim$edge.length * 2

# Node area table (only bifurcating nodes) ----
node_area <- data.frame(
  area = c("ABC", "B", "C", "ABC"),
  row.names = paste0("N", 6:9)
)

# Identify insertion points using ggtree ----
gdata <- ggtree::ggtree(tree_sim)$data
ins_data <- gdata %>% dplyr::filter(node %in% c(1, 3, 7, 8))

# Define inserted nodes with areas ----
inserts <- tibble::tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.5, 0.5, 1, 0.5),
  node_area = c("D", "BC", "AB", "BC")
)

# Insert nodes ----
result <- insert_nodes(tree_sim, inserts, node_area = node_area)
tree_out <- result$phylo
node_area_out <- result$node_area

# Tree visualization ----
tree_viz <- tree_out
tree_viz$node.label <- node_area_out$area
plot(tree_viz, show.node.label = TRUE)
ape::nodelabels(adj = 1.2)
ape::axisPhylo()

# Community presence/absence matrix ----
W <- matrix(
  c(
    0,0,0,1,1,
    1,1,1,1,1,
    0,0,1,1,1,
    1,1,1,0,1,
    1,0,0,1,1
  ),
  nrow = 5, ncol = 5, byrow = TRUE,
  dimnames = list(
    c("Comm 1", "Comm 2", "Comm 3", "Comm 4", "Comm 5"),
    paste("s", 1:5, sep = "")
  )
)

# Biogeographic regions ----
biogeo <- data.frame(Ecoregion = c("A", "B", "B", "C", "D"))




# TESTS -------------------------------------------------------------


test_that("calc_insitu_diversification handles invalid type", {
  expect_error(
    calc_insitu_diversification(
      W, tree_out, node_area_out, biogeo, type = "wrong"
    ),
    '"type" must be "equal.splits" or "fair.proportion"'
  )
})

test_that("calc_insitu_diversification returns expected structure", {
  res <- calc_insitu_diversification(
    W, tree_out, node_area_out, biogeo, type = "equal.splits"
  )
  
  # Must be a list with specific elements
  expect_type(res, "list")
  expect_true(all(c("jetz_site_sp", "jetz_comm_mean",
                    "insitu_site_sp", "insitu_comm_mean",
                    "prop_site_sp", "prop_comm_mean") %in% names(res)))
  
  # Dimensions must match W
  expect_equal(dim(res$jetz_site_sp), dim(W))
  expect_equal(dim(res$insitu_site_sp), dim(W))
  expect_equal(dim(res$prop_site_sp), dim(W))
  expect_equal(length(res$jetz_comm_mean), nrow(W))
  expect_equal(length(res$insitu_comm_mean), nrow(W))
  expect_equal(length(res$prop_comm_mean), nrow(W))
})

test_that("calc_insitu_diversification numerical consistency", {
  res <- calc_insitu_diversification(
    W, tree_out, node_area_out, biogeo, type = "equal.splits"
  )
  
  # prop_comm_mean must equal rowSums(prop_site_sp)/rowSums(jetz_site_sp)
  expect_equal(
    res$prop_comm_mean,
    rowSums(res$prop_site_sp, na.rm = TRUE) /
      rowSums(res$jetz_site_sp, na.rm = TRUE)
  )
  
  # All jetz_comm_mean must be finite
  expect_true(all(is.finite(res$jetz_comm_mean)))
})

test_that("insitu_comm_mean and jetz_comm_mean behave consistently", {
  res <- calc_insitu_diversification(
    W, tree_out, node_area_out, biogeo, type = "equal.splits"
  )
  
  # insitu_comm_mean should be NA or non-negative
  expect_true(all(res$insitu_comm_mean >= 0 | is.na(res$insitu_comm_mean)))
  
  # insitu_comm_mean should equal harmonic mean of insitu_site_sp per row
  expected_insitu <- apply(res$insitu_site_sp, 1, function(x) {
    calc_harmonic_mean(x, na.rm = TRUE, ignore.zero = TRUE)
  })
  expect_equal(res$insitu_comm_mean, expected_insitu)
  
  # jetz_comm_mean should equal harmonic mean of jetz_site_sp per row
  expected_jetz <- apply(res$jetz_site_sp, 1, function(x) {
    calc_harmonic_mean(x, na.rm = TRUE)
  })
  expect_equal(res$jetz_comm_mean, expected_jetz)
})


test_that("prop_comm_mean is computed correctly", {
  res <- calc_insitu_diversification(
    W = W,
    tree = tree_out,
    ancestral.area = node_area_out,
    biogeo = biogeo,
    type = "equal.splits"
  )
  
  # recompute from the matrices
  recomputed_prop <- rowSums(res$prop_site_sp, na.rm = TRUE) /
    rowSums(res$jetz_site_sp, na.rm = TRUE)
  
  expect_equal(res$prop_comm_mean, recomputed_prop, tolerance = 1e-8)
})

