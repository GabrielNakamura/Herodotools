
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


function(x,
         phy, 
         sp.bin = "Sturges", # eu mudaria para breaks aqui 
         marginal = FALSE, 
         lik.threshold = FALSE, 
         threshold = 0.7, 
         compute.fields = F, 
         plot.results = FALSE){
  if(any(x > 1) == TRUE){
    x <- ifelse(x >= 1, 1, 0)
  }
  
  # Enter and organize data:
  match <- picante::match.phylo.comm(phy, x)
  x <- match$comm
  phy <- ape::makeNodeLabel(phangorn::nnls.tree(cophenetic(match$phy),
                                                match$phy,rooted=TRUE), prefix = "Node")
  root.age <- max(cophenetic(phy))/2
  
  # Extract species by nodes matrix
  spp.nodes <- suppressWarnings(t(spp_nodes(tree = phy)))
  
  
  # Compute node ages:
  ages <- round(as.matrix(root.age - ape::node.depth.edgelength(phy)), 5)
  ages[1:length(phy$tip.label),] <- 0.00001 # ages for tip 
  colnames(ages) <- "age"
  rownames(ages) <- c(phy$tip.label, colnames(spp.nodes))
  
  # Compute "LTT":
  h <- hist(ages, breaks = sp.bin, plot = F)
  breaks <- h$breaks[-1]
  n.breaks <- length(breaks)
  mid.time <- as.vector(h$mids)
  
  # Compute weights for species and internal nodes given "LTT":
  weights <- 1/h$counts
  isna.weights <- which(is.na(ifelse(weights == Inf, NA, weights)))
  if(length(isna.weights) == 0){
    weights <- weights
  } else{
    weights <- weights[-isna.weights]
  }
  spp.weight.class <- list()
  spp.weight.class[[1]] <- which(ages <= breaks[1])
  for(b in 1:n.breaks){
    spp.weight.class[[b+1]] <- which(ages > breaks[b] & ages <= breaks[b + 1])
  }
  spp.weight.class <- spp.weight.class[-length(spp.weight.class)]
  spp.weight.class[lengths(spp.weight.class) == 0] <- NA
  isna.spp.weight.class<-which(is.na(spp.weight.class))
  if(length(isna.spp.weight.class) == 0){
    spp.weight.class.clean <- spp.weight.class
  } else{
    spp.weight.class.clean <- rlist::list.remove(spp.weight.class,
                                                 range = which(is.na(spp.weight.class)))
  }
  
  weight.mat <- matrix(NA, 
                       nrow = nrow(ages), 
                       ncol = length(weights), 
                       dimnames = list(rownames(ages), round(weights,3)
                       )
  )
  for(w in 1:nrow(weight.mat)){
    # w = 1
    for(y in 1:length(weights)){
      # y = 1
      weight.mat[w, y] <- round(ifelse(w%in%spp.weight.class.clean[[y]], yes = weights[y], no = 0), 3)
    }
  }
  weight.vec <- vector(length = nrow(weight.mat))
  for(w in 1:nrow(weight.mat)){
    weight.vec[w] <- sum(weight.mat[w, ])
  }
  
  # Run Ancestral Area Reconstruction:
  node.list <- vector(mode = "list", length = nrow(x))
  for(i in 1:nrow(x)){
    # i = 15
    node.list[[i]] <- ape::ace(t(x)[,i], phy, scaled = F, type = "discrete", marginal = marginal)$lik.anc
  }
  threshold <- threshold - (1 - threshold) 
  if(lik.threshold==TRUE){
    node.list.lik<-node.list
    for (k in 1:length(node.list)){
      for (l in 1:phy$Nnode){
        if(abs(node.list[[k]][l, 1] - node.list[[k]][l, 2]) < threshold){
          node.list.lik[[k]][l,] <- 0
        } else {
          node.list.lik[[k]][l,] <- node.list[[k]][l,]
        }
      }
    }
  } else {
    node.list.lik<-node.list
  }
  node.anc.area <- matrix(NA, 
                          nrow = phy$Nnode, 
                          ncol = nrow(x),
                          dimnames = list(phy$node.label, 
                                          rownames(x)
                          )
  )
  for (j in 1:length(node.list.lik)){
    node.anc.area[,j] <- apply(node.list.lik[[j]], 1, which.max)
  }
  node.anc.area <- ifelse(node.anc.area == 1, 0, 1) 
  
  # Define joint occurrence matrix:
  joint.x <- as.matrix(cbind(ifelse(x >= 1, 1, 0), t(node.anc.area)))
  
  # Compute node age by sites matrix:
  age.anc.area <- matrix(0, 
                         nrow(joint.x),
                         ncol(joint.x),
                         dimnames = list(rownames(joint.x), 
                                         colnames(joint.x))
  )
  for (p in 1:nrow(joint.x)){
    age.anc.area[p,] <- ages*t(joint.x)[,p]
  }
  
  # Compute results:
  res.dens <- matrix(NA,
                     nrow = nrow(age.anc.area), 
                     ncol = 7, 
                     dimnames = list(rownames(age.anc.area),c("N_species","N_ancestral_nodes",
                                                              "Density.Peak.Myr","Skewness","Lowest.Distance.To.Peak",
                                                              "Highest.Distance.To.Peak","Peak.Range")
                     )
  )
  
  res.dens[, 1] <- rowSums(x) # number of present day species
  res.dens[, 2] <- rowSums(joint.x[ , (length(phy$tip.label) + 1):ncol(joint.x)]) # Number of ancestors at each site
  for(l in 1:nrow(age.anc.area)){
    #l = 4
    nz.node.ages <- which(age.anc.area[l, ] > 0)
    weights.nz <- weight.vec[nz.node.ages]/sum(weight.vec[nz.node.ages])
    age.anc.area.nz <- age.anc.area[l, nz.node.ages]
    if(length(age.anc.area.nz) >= 2){
      dens.area <- density(age.anc.area.nz, from = 0, to = root.age, weights = weights.nz)
      res.dens[l,3] <- dens.area$x[which(dens.area$y == max(dens.area$y))] # age with most number of species
      res.dens[l,4] <- ifelse(res.dens[l,2] == 0, NA, moments::skewness(dens.area$y))  
      hdi.dens <- suppressWarnings(HDInterval::hdi(dens.area, allowSplit = F, 0.9))
      res.dens[l,5] <- as.numeric(hdi.dens[1]) # lowest distance to peak (myr)
      res.dens[l,6] <- as.numeric(hdi.dens[2]) # highest distance to peak (myr)
      res.dens[l,7]<-res.dens[l,6]-res.dens[l,5] # peak range
    } else{
      res.dens[l,3:7] <- NA
    }
  } 
  
  # Compute site number of nodes x time period:
  div.nodes <- matrix(NA, nrow = nrow(joint.x), 
                      ncol = length(spp.weight.class.clean), 
                      dimnames = list(rownames(joint.x), 
                                      c(mid.time[-which(is.na(spp.weight.class))])
                      )
  )
  names_time_slice <- vector(length = (length(breaks) -1))
  for(i in 1:(length(breaks) -1)){
    names_time_slice[i] <- paste(breaks[i], breaks[i+1], sep = "-")
  }
  for(d in 1:length(spp.weight.class.clean)){
    if(length(spp.weight.class.clean[[d]])>1){
      div.nodes[ , d] <- rowSums(joint.x[ , spp.weight.class.clean[[d]]])
    } else{
      div.nodes[, d] <- joint.x[ , spp.weight.class.clean[[d]]]
    }
  }
  
  # Compute species fields:
  if(compute.fields == TRUE){
    hist.fields.div <- list()
    hist.fields.anc.nodes <- list()
    hist.fields.peak <- list()
    for(f in 1:ncol(x)){
      sites <- which(x[ , f]>0)
      hist.fields.div[[f]] <- res.dens[sites, 1]
      hist.fields.anc.nodes[[f]]<-res.dens[sites,2]
      hist.fields.peak[[f]]<-res.dens[sites,3]
    }
    Species.Fields<-list(Diversity.Fields=hist.fields.div,
                         Ancient.Node.Diversity.Fields=hist.fields.anc.nodes,
                         Diversity.Peak.Time.Fields=hist.fields.peak)
  } else{
    Species.Fields<-"Tell me to compute them, dude!!!"
  }   
  
  
  # Define outputs: 
  Res<-list(Phylogeny=phy,
            Root.Age=root.age,
            Species.per.node.matrix=spp.nodes,
            Node.ages=ages,
            Per.node.ancestral.area=node.anc.area,
            Diversity.Through.Time=div.nodes,
            Cell.Metrics=res.dens,
            Species.Fields)
  return(Res)
  
  # Spatial objects
  if(plot.results == TRUE){
    if(length(coords) == 0){
      stop("Provide a matrix with coordinates to plot the results")
    }
    if(dim(coords) > 2){
      stop("The coordinate matrix should contain two columns (lat and long)")
    }
    shp_earth <- rnaturalearth::ne_countries(returnclass = "sf")
    box <- c(xmin = min(coords[, 1]), xmax = max(coords[, 1]), ymin = min(coords[, 2]), ymax = max(coords[, 2]))
    shp_data <- sf::st_crop(shp_earth, sf::st_bbox(box))
    cell_metrics_df <- data.frame(Res$Cell.Metrics, ID_comm = rownames(comm))
    sf_cell_metrics <-
      shp_data %>% 
      dplyr::left_join(cell_metrics_df)
  }
  
}


# testing plot ada fiunction ----------------------------------------------

ada.res <- Res
dim(Res$Cell.Metrics)
coords <- coord_tyranidae
dim(coord_tyranidae)
match(rownames(coord_tyranidae), rownames(Res$Cell.Metrics))
resolution <- 1
readRDS(file = "grid_tyranidae.rds")
grid <- gridded

plot_ada <- 
  function(ada.res, grid,  patterns, palette = "SunsetDark"){
    ada.res <- ada.res$Cell.Metrics
    # box <- c(xmin = min(coords[, 1]), xmax = max(coords[, 1]), ymin = min(coords[, 2]), ymax = max(coords[, 2]))
    extend_grid <- raster::extend(grid)
    r <- raster::raster(vals = NA, xmn = extend_grid[1],
                        xmx = extend_grid[2],
                        ymn = extend_grid[3],
                        ymx = extend_grid[4]
    )
    # r <- raster::raster(vals = NA, xmn = -170.2166 , xmx = -13.21288, ymn = -55.37714, ymx = 83.6236,
    #                    resolution = resolution)
    #r <- raster::raster(vals = NA, 
    #                    xmn =  min(coords[, 1]), 
    #                    xmx = max(coords[, 1]), 
    #                    ymn = min(coords[, 2]), 
    #                    ymx = max(coords[, 2]),
    #                    resolution = resolution)
    cell.r <- raster::cellFromXY(r, coords[rownames(ada.res),])
    values_cell <- rep(NA, raster::ncell(r))
    names(values_cell) <- 1:raster::ncell(r)
    val.cells <- 1:raster::ncell(r) %in% cell.r
    values_cell[val.cells] <- ada.res[, 2]
    r.n_nodes <- raster::setValues(r, values = values_cell)
    projcrs <- "+proj=robin"
    projection(r.n_nodes) <- projcrs
    df_r_nodes <- as.data.frame(r.n_nodes, xy = T)
    spatial_plot <- 
      ggplot() +
      geom_raster(data = na.omit(df_r_nodes), aes(x = x, y = y, fill = layer)) +
      rcartocolor::scale_fill_carto_c(palette = palette
      )
    
  }





