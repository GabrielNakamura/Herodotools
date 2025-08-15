data(toy_treeEx)
tree_sim <- toy_treeEx
tree_sim$node.label <- NULL

# Create node area table (only bifurcating nodes)
node_area <- data.frame(
  area = c("ABC", "B", "C", "ABC"),
  row.names = paste0("N", 6:9)
)

# create insertion of nodes to simulate for change in area in the branch 
# Use ggtree to identify insertion points 

gdata <- ggtree::ggtree(tree_sim)$data
ins_data <- gdata %>% dplyr::filter(node %in% c(7, 8))
ins_data <- tibble::add_row(ins_data, ins_data[1, ]) %>% dplyr::arrange(node)

# Define inserted nodes with areas
inserts <- tibble::tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.25, 0.75, 0.25),
  node_area = c("AB", "B", "BC")
)

# Insert nodes 
result <- insert_nodes(tree_sim, inserts, node_area = node_area)
tree_out <- result$phylo
node_area_out <- result$node_area

# tree visualization
tree_viz <- tree_out
tree_viz$node.label <- node_area_out$area

plot(tree_viz, show.node.label = T)
ape::nodelabels(adj = 1.2)


W <- matrix(
  c(
    0,1,0,1,1,
    1,1,1,1,1,
    1,1,0,1,1,
    0,1,0,1,1),
  nrow= 4, ncol= 5, byrow = T,
  dimnames=list(c("Comm 1", "Comm 2", "Comm 3", "Comm 4"),c(paste("s", 1:5, sep=""))))


biogeo <- data.frame(Ecoregion= c("A", "B", "C", "D"))

# comparando os resultados as funcoes divergem para a DR insitu
# investigando os resultados descobir que a funcao antiga pode estar calculando errado o ED
# ela está incluindo o ramo do ancestral na conta, qndo deveria contar apenas a parte descendente


insitu2 <- calc_insitu_diversification2(W = W,
                            tree = tree_sim, 
                            ancestral.area = node_area,
                            biogeo = biogeo, 
                            diversification = c("jetz"))

calc_insitu_diversification(W = W,
                            tree = tree_sim, 
                            ancestral.area = node_area,
                            biogeo = biogeo, 
                            diversification = c("jetz"))



calc_insitu_diversification2(W = W,
                             tree = tree_sim, 
                             ancestral.area = node_area,
                             biogeo = biogeo, 
                             diversification = c("freckleton"))


# testes sobre media hamonica para DR in situ
harmonic_mean <- function(x) {
  x <- x[!(x==0 | x==0.00001)]

  length(x) / sum(1 / x)
}

harmonic_mean(insitu2$Jetz_species_site[2,])

harmonic_mean(insitu2$insitu_Jetz_species_sites[2,])


calc_node_density(tree = tree_sim, 
                  ancestral.area = node_area,
                  current.area = "D")

