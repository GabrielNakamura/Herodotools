# tests/testthat/test-calc_bsm.R

test_that("calc_bsm runs and returns expected output structure", {
  # Load example data from the Herodotools package
  data("akodon_sites")
  data("akodon_newick")
  
  site_xy <- akodon_sites %>%
    dplyr::select(LONG, LAT)
  
  akodon_pa <- akodon_sites %>%
    dplyr::select(-LONG, -LAT)
  
  spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
  akodon_pa_tree <- akodon_pa[, spp_in_tree]
  
  # Load ancestral reconstruction
  load(file = system.file("extdata", "resDEC_akodon.RData", package = "Herodotools"))
  
  n_maps <- 3
  
  # Run BSM with small number of mappings for test speed
  bsm <- calc_bsm(
    BioGeoBEARS.data = resDEC,
    phyllip.file = system.file("extdata", "geo_area_akodon.data", package = "Herodotools"),
    tree.path = system.file("extdata", "akodon.new", package = "Herodotools"),
    max.maps = 10,
    n.maps.goal = 3,
    seed = 123
  )
  
  # Basic structure checks
  expect_type(bsm, "list")
  expect_length(bsm, 2)
  expect_named(bsm)
  expect_true(
    all(c("RES_clado_events_tables", "RES_ana_events_tables" ) %in% names(bsm))
    )
  expect_true(length(bsm$RES_clado_events_tables) >= n_maps)
  
  # Check that the number of stochastic maps returned matches n.maps.goal 
  expect_true(length(bsm$RES_clado_events_tables) >= n_maps)
  expect_true(length(bsm$RES_ana_events_tables) >= n_maps)
  
})








