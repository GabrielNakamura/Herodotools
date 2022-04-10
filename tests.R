
# general installation ----------------------------------------------------

devtools::install_github("GabrielNakamura/Rrodotus", ref = "main", force = TRUE)
library(Rrodotus)
data("comm_data")
data("biogeo")
data("node_biogeo")
data("tree_aves")



# testing DB metrics function ---------------------------------------------

test_leandro <- DivB_metrics(W = comm_data,
                             tree = tree_aves,
                             ancestral.area = node_biogeo,
                             biogeo = biogeo,
                             diversification = c("jetz", "freck"),
                             PD = TRUE,
                             PE = TRUE,
                             age.arrival = TRUE,
                             age.no.ancestor = NA,# 'half.edge' or numeric()
                             dispersal.from = TRUE,
                             ED.type = "equal.splits"
)



# test ada function -------------------------------------------------------

# akodon data

comm_akodon <- read.table("comm_akodon.txt", header = TRUE)
coords_akodon <- read.table("coord_akodon.txt", header = TRUE)
size_akodon <- read.table("size_akodon.txt", header = TRUE)
phy_akodon <- ape::read.nexus("tree_akodon.nexus")


# simulating communities

nsp <- 100
ncomm <- 20
comm <- matrix(rpois(nsp*ncomm, 1), nrow = ncomm, ncol = nsp,
               dimnames = list(paste("comm", 1:ncomm, sep = "_"),
                               paste("s", 1:nsp, sep = ""))
)
phy <- geiger::sim.bdtree(b = 1, d = 0, n = nsp)
x <- comm 
phy <- phy
sp.bin = "Sturges"
marginal = FALSE
lik.threshold = FALSE
threshold = 0.7
compute.fields = F

# akodon communities
phy <- phy_akodon
x <- comm_akodon 
sp.bin = "Sturges"
marginal = FALSE
lik.threshold = FALSE
threshold = 0.7
compute.fields = F
plot.results = F
coords = coords_akodon

# tyranidae communities
phy <- phylo_tyranidae
x <- comm_tyranidae 
sp.bin = "Sturges"
marginal = FALSE
lik.threshold = FALSE
threshold = 0.7
compute.fields = F
plot.results = F
coords = coord_tyranidae
shp_tyranidae # spatial polygons for Tyranidae

test_ada_tyranidae <- ada(x = comm_tyranidae,
                          phy = phylo_tyranidae, 
                          sp.bin = "Sturges",
                          marginal = FALSE, 
                          lik.threshold = FALSE, 
                          threshold = 0.7, 
                          compute.fields = FALSE)
saveRDS(object = test_ada_tyranidae, "res_ada_tyra.rds")

library(magrittr)
library(ggplot2)
library(sf)
library(raster)

res_plot_test <- plot_ada(ada.res = res_ada, 
                          grid = shp_tyranidae, 
                          coords = coord_tyranidae, 
                          resolution = 1, 
                          patterns = "all",
                          color_palette = "SunsetDark")
dev.new()
res_plot_test$Nnodes
res_plot_test$HighDistPeak
dev.new()

