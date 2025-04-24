get_path_between_nodes <- function(gdata, start, end) {
  path <- c(start)
  current <- start
  
  while (current != end) {
    parent <- gdata$parent[gdata$node == current]
    if (is.na(parent)) {
      stop("End node is not an ancestor of the start node.")
    }
    path <- c(path, parent)
    current <- parent
  }
  
  final_path <- path[-length(path)]
  # Return rows in gdata for the nodes along the path
  gdata %>% filter(node %in% final_path)
}



# gdata <- ggtree(tree_out)$data
# 
# # Get path from tip 3 to node 6
# path_df <- get_path_between_nodes(gdata, start = 3, end = 7)
# print(path_df)
# 
# inserts
# 
# round(bl_orig, 10) == round(sum(path_df$branch.length), 10)
