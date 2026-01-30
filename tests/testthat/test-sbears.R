
# generating data for tests -----------------------------------------------

# phylogenetic tree
phylo <- geiger::sim.bdtree(n = 10, seed = 42)
phylo <- ape::makeNodeLabel(phy = phylo)

# community composition matrix
comm <- 
  matrix(sample(c(0, 1), size = 10*20, replace = TRUE),
         nrow = 20, 
         ncol = 10, 
         dimnames = list(paste("comm", 1:20), phylo$tip.label)
  )

# coordinates - this is necessary for "disperal_assembly" algorithm
xy_coords <- 
  matrix(runif(1:10), 
         nrow = nrow(comm),
         ncol = 2, 
         dimnames = list(rownames(comm), c("lng", "lat")
         )
  )

# running sbears with "single_site" algorithm
out_sbears_single_node_comp <- 
  calc_sbears(x = comm,
              phy = phylo, 
              coords = xy_coords, 
              method = "single_site", 
              compute.node.by.sites = TRUE)

# matrix containing ancestral reconsturction, nodes are rows and columns are sites
anc_reconstruction_single <- out_sbears_single_node_comp$reconstruction
node_comp_single <- out_sbears_single_node_comp$site_node_composition

# running sbears with "dispersal_assembly" algorithm
out_sbears_disp_node_comp <- 
  calc_sbears(x = comm, 
              phy = phylo, 
              coords = xy_coords, 
              method = "disp_assembly", 
              w_slope = 5, 
              min_disp_prob = 0.8, 
              compute.node.by.sites = TRUE)

anc_reconstruction_disp <- out_sbears_disp_node_comp$reconstruction
node_comp_disp <- out_sbears_disp_node_comp$site_node_composition

# testing output format and reliability -----------------------------------

testthat::test_that("difference in output between single site and dispersal", {
  expect_equal(node_comp_single, node_comp_disp) 
  expect_equal(dim(anc_reconstruction_single), dim(anc_reconstruction_disp))
  expect_true(sum(anc_reconstruction_disp != anc_reconstruction_single) >= 1, TRUE)
})


testthat::test_that("output is in the right format and dimensions", {
  
  # checking if the dimensions between outputs are the same for both algorithms
  expect_equal(all.equal(dim(t(anc_reconstruction_single)), 
                         dim(node_comp_single)), 
               TRUE)
  expect_equal(all.equal(dim(t(anc_reconstruction_disp)), 
                         dim(node_comp_disp)), 
               TRUE)
  
  # checking if dimensions of output with different algorithms but same data are the same
  expect_equal(all.equal(dim(anc_reconstruction_disp), 
                         dim(anc_reconstruction_single)), 
               TRUE)
  
  # checking if the number of communities are the same in the input and outputs
  expect_equal(nrow(comm) == ncol(anc_reconstruction_disp) &
                 nrow(comm) == ncol(anc_reconstruction_single),
               TRUE)
  expect_equal(nrow(node_comp_disp) == nrow(comm) & 
                 nrow(node_comp_single == nrow(comm)), 
               TRUE)
  
  # checking if all the nodes of the tree are in output
  expect_equal(nrow(anc_reconstruction_disp) == phylo$Nnode & 
                 nrow(anc_reconstruction_single) == phylo$Nnode, TRUE)
  
  # check if the function is returning probabilities 
  expect_true(all.equal(anc_reconstruction_disp <= 1, 
                            anc_reconstruction_single <= 1),
                  TRUE)
  
  # checking if the sites names are preserved in output 
  expect_true(all.equal(colnames(anc_reconstruction_disp) == rownames(comm), 
                colnames(anc_reconstruction_single) == rownames(comm)), TRUE)
  
})
