
# arvore para exemplo
tree_topo <- "(((s1:1,s2:1):1,(s3:1,s4:1):1):5,(((s5:1,s6:1):1,(s7:1,s8:1):1):1,s9:3):4);"
tree_test <- ape::read.tree(text = tree_topo)
tree_test <- ape::makeNodeLabel(tree_test, prefix = "N", method = "number")

# resultado da reconstrucao 
comm_obs <- matrix(1, nrow = 1, ncol = 4, dimnames = list("comm_1", c("s1", "s2", "s8", "s9")))
nodes_reconstruction <- data.frame(comm = "comm_1", nodes = c("N3", "N2", "N1"))
nodes_comm <- data.frame(comm = "comm_1", nodes = c("N3", "N1", "N5"))

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

# PD
table_tree_potential <- tidytree::as_tibble(tree_potential)
library(dplyr)
library(tidyr)
table_tree_potential2 <- 
  table_tree_potential %>% 
  mutate(pres = ifelse(label %in% colnames(comm_obs), "pres", "abs")) %>% 
  mutate(ancestor = table_tree_potential$label[table_tree_potential$parent]) %>% 
  mutate(descendant = table_tree_potential$label[table_tree_potential$node]) %>% 
  mutate(partition1 = gsub(pattern = ".*_", replacement = "", x = ancestor)) %>% 
  mutate(partition2 = gsub(pattern = ".*_", replacement = "", x = descendant)) %>% 
  mutate(partition.IS = ifelse(partition1 == "IS" | partition2 == "IS", "IS", NA)) %>% 
  mutate(partition.IM = ifelse(partition1 == "IM" & pres == "pres" | partition1 == "ESD" & pres == "pres", "IM", NA)) %>% 
  mutate(partition.EM = ifelse(partition1 == "EM" & partition2 == "ESD", "EM", NA)) %>% 
  mutate(partition.ESD = ifelse(is.na(partition.IS) == TRUE & is.na(partition.IM) == TRUE & is.na(partition.EM) == TRUE, "ESD", NA))
  
# calculo PD

PDinsitu <- 
  table_tree_potential2 %>% 
  filter(partition.IS == "IS") %>% 
  select(branch.length) %>% 
  sum(na.rm = T)
PDimmigration <- 
  table_tree_potential2 %>% 
  filter(partition.IM == "IM") %>% 
  select(branch.length) %>% 
  sum(na.rm = T)

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

PDtotal <- PDinsitu + PDimmigration + PDemmigration + PDexsitu

data_res <- 
  data.frame(partition = c("PDinsitu", "PDimmigration", "PDemmigration", "PDexsitu", "PDtotal"), 
             value = c(PDinsitu, PDimmigration, PDemmigration, PDexsitu, PDtotal), community = comm_names[i])
