
# arvore para exemplo
tree_topo <- "(((s1:1,s2:1):1,(s3:1,s4:1):1):5,(((s5:1,s6:1):1,(s7:1,s8:1):1):1,s9:3):4);"
tree_test <- ape::read.tree(text = tree_topo)
tree_test <- ape::makeNodeLabel(tree_test, prefix = "N", method = "number")

# resultado da reconstrucao 
comm_obs <- matrix(1, nrow = 1, ncol = 5, dimnames = list("comm_1", c("s1", "s2", "s5", "s8", "s9")))
nodes_reconstruction <- data.frame(comm = "comm_1", nodes = c("N3", "N2", "N7", "N8", "N5"))
nodes_comm <- data.frame(comm = "comm_1", nodes = c("N3", "N1", "N6", "N5"))

# encontrando a arvore potencial
node_sequence <- gsub(pattern = "N", replacement = "", unique(c(nodes_reconstruction$nodes, nodes_comm$nodes)))
data_nodes <- data.frame(nodes = unique(c(nodes_reconstruction$nodes, nodes_comm$nodes)), node_sequence = as.numeric(node_sequence))
basal_pos <- which(data_nodes$node_sequence == min(data_nodes$node_sequence))
tree_potential <- ape::extract.clade(phy = tree_test, node = data_nodes[basal_pos, 1])

# ISDiv
nodes_insitu <- intersect(nodes_reconstruction$nodes, nodes_comm$nodes)
tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(nodes_insitu, tree_potential$node.label)], "IS", sep = "_")

# IM
node_immigration <- setdiff(nodes_comm$nodes, nodes_reconstruction$nodes)
tree_potential$node.label[match(node_immigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_immigration, tree_potential$node.label)], "IM", sep = "_")

# EM
node_emmigration <- setdiff(nodes_reconstruction$nodes, nodes_comm$nodes)
tree_potential$node.label[match(node_emmigration, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_emmigration, tree_potential$node.label)], "EM", sep = "_")

# ESDiv
node_exsitu <- tree_potential$node.label[-grep(pattern = "_", tree_potential$node.label)] 
tree_potential$node.label[match(node_exsitu, tree_potential$node.label)] <- paste(tree_potential$node.label[match(node_exsitu, tree_potential$node.label)], "ESD", sep = "_")

# table node combinations
table_tree_potential <- tidytree::as_tibble(tree_potential)
library(dplyr)
library(tidyr)

table_tree_potential2 <- 
  table_tree_potential %>% 
  mutate(pres = ifelse(label %in% colnames(comm_obs), "pres", "abs")) %>% 
  mutate(ancestor = table_tree_potential$label[table_tree_potential$parent]) %>% 
  mutate(descendant = table_tree_potential$label[table_tree_potential$node]) %>% 
  mutate(ancestor1 = gsub(pattern = ".*_", replacement = "", x = ancestor)) %>% 
  mutate(descendent1 = gsub(pattern = ".*_", replacement = "", x = descendant)) %>% 
  mutate(partition.IS = ifelse(ancestor1 == "IS" & descendent1 == "IS" |
                                 ancestor1 == "IS" & descendent1 == "EM" |
                                 ancestor1 == "EM" & descendent1 == "IS" | 
                                 ancestor1 == "EM" & descendent1 == "EM" |
                                 ancestor1 == "IS" & pres == "pres" |
                                 ancestor1 == "EM" & pres == "pres",
                               "IS", NA)) %>% 
  mutate(partition.IM = ifelse(ancestor1 == "IM" & descendent1 == "IS" |
                                 ancestor1 == "IM" & descendent1 == "EM",
                               "IM", NA)) %>% 
  mutate(partition.EM = ifelse(ancestor1 == "IS" & descendent1 == "IM" |
                                 ancestor1 == "EM" & descendent1 == "ESD" |
                                 ancestor1 == "EM" & descendent1 == "IM" | 
                                 ancestor1 == "EM" & pres == "abs" & descendent1 != "IS", 
                               "EM", NA)) %>% 
  mutate(partition.ESD = ifelse(ancestor1 == "ESD" & descendent1 == "ESD" |
                                  ancestor1 == "ESD" & pres == "abs", "ESD", NA)) %>% 
  mutate(partition.unknown = ifelse(ancestor1 == "IM" & descendent1 == "IM", "unknown", NA))




# PD decomposition

PDinsitu <- 
  table_tree_potential2 %>% 
  filter(partition.IS == "IS") %>% 
  select(branch.length) %>% 
  sum(na.rm = T)
PDimmigration <- 
  table_tree_potential2 %>% 
  filter(partition.IM == "IM") %>% 
  select(branch.length) %>% 
  mutate(qnt0 = lapply(branch.length, function(x) quantile(seq(from = 0.0001, to = x, by = 0.001)))[1]) %>% 
  mutate(qnt25 = lapply(branch.length, function(x) quantile(seq(from = 0.0001, to = x, by = 0.001)))[[1]][2]) %>% 
  mutate(qnt50 = lapply(branch.length, function(x) quantile(seq(from = 0.0001, to = x, by = 0.001)))[[1]][3])

  sum(na.rm = T)

unlist(lapply(5, function(x) quantile(seq(from = 0.0001, to = x, by = 0.001))))
data("iris")
iris %>% 
  mutate(Max.Len= purrr::pmap_dbl(list(Sepal.Length, Petal.Length), max))

table_tree_potential2 %>% 
  filter(partition.IM == "IM") %>% 
  select(branch.length) %>% 
  mutate(qnt0 = purrr::pmap_dbl(list(branch.length), quantile(seq(from = 0.0001,  by = 0.001))))

PDemmigration <- 
  table_tree_potential2 %>% 
  filter(partition.EM == "EM") %>% 
  select(branch.length) %>% 
  sum(na.rm = TRUE)
PDexsitu <- 
  table_tree_potential2 %>% 
  filter(partition.ESD == "ESD") %>% 
  select(branch.length) %>% 
  sum(na.rm = TRUE)

PDunknown <- 
  table_tree_potential2 %>% 
  filter(partition.unknown == "unknown") %>% 
  select(branch.length) %>% 
  sum(na.rm = TRUE)


PDtotal <- PDinsitu + PDimmigration + PDemmigration + PDexsitu + PDunknown

data_res <- 
  data.frame(partition = c("PDinsitu", "PDimmigration", "PDemmigration", "PDexsitu", "PDtotal", "PDunknown"), 
             value = c(PDinsitu, PDimmigration, PDemmigration, PDexsitu, PDtotal, PDunknown), community = "comm_1")


PDinsitu + PDimmigration

picante::pd(matrix(1, nrow = 1, ncol = 5, dimnames = list("comm_1", c("s1", "s2", "s5", "s8", "s9"))), tree = tree_potential)
