
# generating data for test ------------------------------------------------

# phylogenetic tree
 set.seed(42)
 phylo <- ape::rcoal(n = 10)
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
 out_sbears <- 
  calc_sbears(x = comm,
              phy = phylo, 
              coords = xy_coords, 
              method = "single_site")
              
 # sbears object can be used directly in PD_decomposition function
 
 out_decomp <- 
   PD_decomposition(comm = comm,
                    sbears.obj = out_sbears, 
                    phy = phylo, 
                    threshold = 0.5)
 
 # Data frame with results from PD decomposition
 out_pd_decomposition_potential <- out_decomp$decomp_potential$PD_decomposition

 out_pd_decomposition_obs <- out_decomp$decomp_faith$PD_decomposition
 
 PD_faith_total <- 
   out_pd_decomposition_obs |> 
   dplyr::group_by(community) |> 
   filter(partition == "PDtotal") |>
   pull(value)

PD_faith_insitu <- 
  out_pd_decomposition_obs |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDinsitu") |>
  pull(value)

PD_faith_exsitu <- 
  out_pd_decomposition_obs |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDexsitu") |>
  pull(value)

PD_faith_emigration <- 
  out_pd_decomposition_obs |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDemigration") |>
  pull(value)

PD_faith_immigration <- 
  out_pd_decomposition_obs |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDimmigration") |>
  pull(value)

# components for estimated PD from sbears

PD_faith_total_potential <- 
  out_pd_decomposition_potential |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDtotal") |>
  pull(value)

PD_faith_insitu_potential <- 
  out_pd_decomposition_potential |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDinsitu") |>
  pull(value)

PD_faith_exsitu_potential <- 
  out_pd_decomposition_potential |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDexsitu") |>
  pull(value)

PD_faith_emigration_potential <- 
  out_pd_decomposition_potential |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDemigration") |>
  pull(value)

PD_faith_immigration_potential <- 
  out_pd_decomposition_potential |> 
  dplyr::group_by(community) |> 
  filter(partition == "PDimmigration") |>
  pull(value)



# testing output format and reliability -----------------------------------

# for observed PD
testthat::test_that("PD total must be always the ceiling", {
  # comparing all communities
  expect_equal(sum(PD_faith_emigration, 
                   PD_faith_exsitu, 
                   PD_faith_immigration, 
                   PD_faith_insitu), sum(PD_faith_total))
  # comparing each community
  expect_equal(colSums(rbind(PD_faith_emigration, 
                   PD_faith_exsitu, 
                   PD_faith_immigration, 
                   PD_faith_insitu)), PD_faith_total)
  
})

# for estimated PD

testthat::test_that("PD potential total must be always the ceiling", {
  # comparing all communities
  expect_equal(sum(PD_faith_emigration_potential, 
                   PD_faith_exsitu_potential, 
                   PD_faith_immigration_potential, 
                   PD_faith_insitu_potential), sum(PD_faith_total_potential))
  # comparing each community
  expect_equal(colSums(rbind(PD_faith_emigration_potential, 
                             PD_faith_exsitu_potential, 
                             PD_faith_immigration_potential, 
                             PD_faith_insitu_potential)), PD_faith_total_potential)
  
})
