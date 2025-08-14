set.seed(4523)
tree_sim <- ape::rcoal(5)

# Create node area table (only bifurcating nodes)
node_area <- data.frame(
  area = c("A", "A", "BC", "AD"),
  row.names = paste0("N", 6:9)
)

# create insertion of nodes to simulate for change in area in the branch 
# Use ggtree to identify insertion points 

gdata <- ggtree::ggtree(tree_sim)$data
ins_data <- gdata %>% dplyr::filter(node %in% c(3, 9, 8))
ins_data <- tibble::add_row(ins_data, ins_data[2, ]) %>% dplyr::arrange(node)

# Define inserted nodes with areas
inserts <- tibble::tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.2, 0.5, 0.9, .3),
  node_area = c("AC", "AB", "BC", "D")
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
    1,0,0,1,0,
    0,0,1,0,1,
    1,0,1,1,1,
    1,1,0,1,1),
  nrow= 4, ncol= 5,
  dimnames=list(c("Comm 1", "Comm 2", "Comm 3", "Comm 4"),c(paste("t", 1:5, sep=""))))


biogeo <- data.frame(Ecoregion= c("A", "A", "C", "D"))

calc_insitu_diversification(W = W,
                            tree = tree_out, 
                            ancestral.area = node_area_out,
                            biogeo = biogeo, 
                            diversification = c("jetz"))

