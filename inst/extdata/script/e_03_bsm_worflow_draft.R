library(devtools)
load_all()
library(ape)
library(BioGeoBEARS)


data("akodon_newick")
data("resDEC")  # BioGeoBEARS model result

# Paths to tree and geography files

tree_path <- system.file("extdata", "akodon.new", package = "Herodotools")
phyllip_path <- system.file("extdata", "geo_area_akodon.data", package = "Herodotools")
max_range_size <- resDEC$inputs$max_range_size


resDEC$inputs$trfn <- tree_path
resDEC$inputs$geogfn <- phyllip_path

## FUNCTION 1 - calc_bsm

bsm_result <- calc_bsm(
   BioGeoBEARS.data = resDEC,
   phyllip.file = phyllip_path,
   tree.path = tree_path,
   max.maps = 50,
   n.maps.goal = 3,
   seed = 1234
 )
 
bsm_result


## new function, not so relevant for the package but good for comparison with
##    biogeobears plots
plot_biogeobears_stochastic_map(resDEC, bsm_result, map_index = 1) 
 

## FUNCTION 2 - get_insert_df

insert_list <- get_insert_df(
   bsm_result,
   phyllip.file = phyllip_path,
   max.range.size = resDEC$inputs$max_range_size
   )

# FUNCTION 3 - get_bsm_node_area

node_area_list <- get_bsm_node_area(
  bsm = bsm_result, 
  BioGeoBEARS.data = resDEC,
  phyllip.file = phyllip_path,
  tree.path = tree_path,
  max.range.size = max_range_size)
  
  

# FUNCTION 4 - insert_nodes
new_nodes_tree <- insert_nodes(
  tree = akodon_newick, 
  inserts = insert_list, 
  node_area = node_area_list)

## All modifications done ---

## VISUALIZATION

plot(new_nodes_tree[[1]]$phylo)
nodelabels()
tiplabels()

new_nodes_tree[[1]]$node_area




library(ggtree)
library(ggplot2)
library(dplyr)
library(tidytree)
library(viridis)


tiprange <- getranges_from_LagrangePHYLIP(lgdata_fn = phyllip_path)

tip_area_df <- data.frame(
  label = rownames(tiprange@df), 
  area = apply(tiprange@df, 1, function(row) {
    paste0(names(row)[which(row == 1)], collapse = "")
  })
) 


new_nodes_tree <- insert_nodes(
  tree = akodon_newick, 
  inserts = insert_list, 
  node_area = node_area_list)





tree <- new_nodes_tree[[1]]$phylo
node_area <- new_nodes_tree[[1]]$node_area

# Start with the base tree using ggtree
p <- ggtree(tree)

# Convert to a data frame to work with positions
tree_data <- as_tibble(p$data)

# Parse node_area to match the ggtree structure
new_tree_data <- node_area %>%
  rownames_to_column("label") %>%
  #mutate(node = as.numeric(gsub("ana_N|N", "", label))) %>%
  add_row(tip_area_df) %>% 
  right_join(tree_data)  


# Choose distinct categorical colors
areas <- sort(unique(na.omit(new_tree_data$area)))
n_areas <- length(areas)

# Better than viridis for categorical: colorspace
biome_colors <- qualitative_hcl(n_areas, palette = "Dark 3")  # or "Set 2", "Dynamic"
names(biome_colors) <- areas

# Plot tree with tip labels
ggplot(new_tree_data) + 
  geom_tree(aes(color = area), size =2) +
  #geom_point(aes(x = x, y = y, color = area), size = 3) +
  geom_point(data = filter(new_tree_data, !grepl("^ana_N", label)),
             aes(x = x, y = y, color = area), size = 3) +
  geom_label(data = filter(new_tree_data, !grepl("^ana_N", label)),
             aes(x = x, y = y, label = area, fill = area), vjust = 0.5, size = 3) +
  scale_color_manual(values = biome_colors) +
  scale_fill_manual(values = biome_colors) +
  geom_tiplab(offset = 1, hjust = 1) +
  
  theme_tree2() +
  ggtitle("Phylogenetic Tree with Added Nodes for Range Change") +
  theme(legend.position = "right")



