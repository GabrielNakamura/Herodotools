
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


res_plot_test <- plot_ada(ada.res = test_ada_tyranidae, 
                          grid = shp_tyranidae, 
                          coords = coord_tyranidae, 
                          resolution = 1, 
                          patterns = "all",
                          color_palette = "SunsetDark")
saveRDS(object = test_ada_tyranidae, "res_ada_tyra.rds")


plot_ada <- 
  function(ada.res,
           grid,
           coords, 
           resolution, 
           patterns = "all",
           color_palette = "SunsetDark")
  {
    
    
    # preparing data and generating raster ------------------------------------
    
    ada.res <- ada.res$Cell.Metrics
    extend_grid <- grid@bbox
    r <- raster::raster(vals = NA, xmn = extend_grid[1, 1],
                        xmx = extend_grid[1, 2],
                        ymn = extend_grid[2, 1],
                        ymx = extend_grid[2, 2], 
                        resolution = resolution
    )
    cell.r <- raster::cellFromXY(r, coords[rownames(ada.res),])
    values_cell <- rep(NA, raster::ncell(r))
    names(values_cell) <- 1:raster::ncell(r)
    val.cells <- 1:raster::ncell(r) %in% cell.r
    
    
    # transforming to sf object -----------------------------------------------
    test_sf <- sf::st_as_sf(raster::rasterToPolygons(r.n_nodes))
    test_sf$ID <- names(values_cell[val.cells])
    test_sf$rich[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 1]
    test_sf$Nnodes[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 2]
    test_sf$PeakDiv[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 3]
    test_sf$Skewness[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 4]
    test_sf$LowDistPeak[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 5]
    test_sf$HighDistPeak[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 6]
    test_sf$PeakRange[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 7]
    
    if(patterns == "all"){
      names_div <- colnames(test_sf)[4:ncol(test_sf)]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            geom_sf(data = test_sf_proj, aes_string(geometry = "geometry", 
                                                    fill = x), 
                    color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    } else{
      names_div <- colnames(test_sf)[4:ncol(test_sf)]
      patterns_plot <- pmatch(patterns, names_div)
      names_div <- names_div[patterns_plot]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            geom_sf(data = test_sf_proj, aes_string(geometry = "geometry", 
                                                    fill = x), 
                    color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    }
    res_plot
  }

