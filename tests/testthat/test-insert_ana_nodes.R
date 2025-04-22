library(ggtree)
library(ape)
library(dplyr)
library(stringr)
library(devtools)
library(patchwork)

load_all()

# create data ----
set.seed(42)
tree <- rcoal(5)

plot(tree)
nodelabels()
tiplabels()


# Create trait table
node_area <- data.frame(
  area = c("A", "A", "BC", "D"), 
  row.names = paste0("N", 6:9)
)


# Find edge to insert on (e.g., tip 3)
gdata <- ggtree(tree)$data
ins_data <- gdata %>% filter(node %in% c(3, 8, 9)) 
ins_data <- add_row(ins_data, ins_data[3,])

# Define insertion
# insert 2 nodes in 2 branches
inserts <- tibble(
  parent = ins_data$parent,
  child = ins_data$node,
  event_time = c(0.02, 0.03, 0.04, .05),
  node_area = c("C", "AB", "AD", "D")
)

# Run function
result <- insert_ana_nodes(tree, inserts, node_area = node_area)

# Now plot
plt_phy <- ggtree(result$phylo) +
  geom_point2(color = "red", size = 3) +
  geom_text2(data = result$data, aes(label = node_area), vjust = -0.3, hjust = -0.3)


plt_dt <- ggtree(result$data) +
  geom_point2(color = "red", size = 3) +
  geom_text2(aes(label = node_area), vjust = -0.3, hjust = -0.3)


plt_phy + plt_dt


# Confirm tree height remains the same
cat("Original tree height:", max(node.depth.edgelength(tree)), "\n")
cat("Modified tree height:", , "\n")



test_that("After node addition, tree has the same total length", {
  expect_equal(
    max(node.depth.edgelength(tree)),
    max(node.depth.edgelength(result$phylo))
    )
})

test_that("Plots from $phylo and $data should produce the same results", {
  expect_equal(plt_phy$data, plt_dt$data)
})
