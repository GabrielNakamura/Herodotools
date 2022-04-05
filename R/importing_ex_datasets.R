
# adding data to run examples ---------------------------------------------

# Akodon communities

usethis::use_data(comm_akodon, Rrodotus)
usethis::use_data(coords_akodon, Rrodotus)
usethis::use_data(size_akodon, Rrodotus)
usethis::use_data(phy_akodon, Rrodotus)

# Tyranidae communities

comm_tyranidae <- read.table(file = "W_harvey.txt", header = TRUE)
coord_tyranidae <- read.table(file = "coords_tyranidae.txt", header = TRUE)
phylo_tyranidae <- ape::read.tree(file = "Tree_TF400Howard_tip_corrected.txt")
shp_tyranidae <- readRDS("grid_tyranidae.rds")
res_ada <- readRDS("res_ada_tyra.rds")

usethis::use_data(comm_tyranidae, Rrodotus)
usethis::use_data(phylo_tyranidae, Rrodotus)
usethis::use_data(coord_tyranidae, Rrodotus)
