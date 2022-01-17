devtools::install_github("GabrielNakamura/Rrodotus", ref = "main", force = TRUE)
library(Rrodotus)
help(package = "Rrodotus") # there are functions with poor documentation
data("comm_data")
data("biogeo")
data("node_biogeo")
data("tree_aves")

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

