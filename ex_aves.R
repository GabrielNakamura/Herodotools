comm_data <- read.table(here::here("dados_exemplo", "community.txt"))
biogeo <- read.table(here::here("dados_exemplo", "biogeo.txt"))
biogeo <- chartr("123456", "ABCDEF", biogeo[, 1])
biogeo <- data.frame(biogeo = biogeo, row.names = rownames(biogeo))
tree_aves <- ape::read.tree(file = here::here("dados_exemplo", "bird_orders.txt"))
node_biogeo <- read.table(here::here("dados_exemplo", "node_biogeo.txt"))


test_leandro <- DivB_metrics(W = comm_data,
         tree = tree,
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
