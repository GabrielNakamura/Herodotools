
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



# Processing Biogeobears data --------------------------------------------------------


load(file = "resBAYAREALIKE.RData")
resBAYAREALIKE
install.packages("rexpokit")
install.packages("cladoRcpp")

library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")
resBAYAREALIKE

# reading data and libraries ----------------------------------------------
load(here::here("output", "resBAYAREALIKE.RData"))


# reading data and libraries ----------------------------------------------
load(here::here("output", "resBAYAREALIKE.RData"))
geogfn = np(paste(here::here("Phyllip_GeoAreas.txt", sep=""))) # Ask renan what is this file
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr <- ape::read.tree(here::here("Tree_TF400Howard_Pruned.tre"))

trtable = BioGeoBEARS::prt(tr, printflag=FALSE)

# extracting probabilites of states/range at each node --------------------
prob_state <- resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_bottom_below_node


# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 3
areas = getareas_from_tipranges_object(tipranges)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = cladoRcpp::rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE) # list with range area for each node and species 

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based)) {
  # i= 175
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}


# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(prob_state)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list

# Write the table to a tab-delimited text file (for Excel etc.)
write.table(range_probabilities, file= here::here("data", "processed", 'range_probabilities.csv'), sep=';', dec=',', row.names=FALSE)

#############################################################################################################
#################################Fazer matrix eco Nodes######################################################
#############################################################################################################

# matrix Node x Biome -----------------------------------------------------

teste_max <- apply(range_probabilities,
                   MARGIN = 1, 
                   function(i) colnames(range_probabilities)[which(i == max(i))
                   ]
) # taking the max probability value

biome <- data.frame(node = 1:length(unlist(teste_max)), biome = unlist(teste_max))

# saving table ------------------------------------------------------------
write.table(biome, file = here::here("data", "processed", 'EcoNodes_harvey.csv'), sep=';', dec=',', row.names=FALSE)

write.table(biome, file = here::here("data", "processed", "Econodes_harvey.txt"))
