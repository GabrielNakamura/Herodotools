devtools::load_all()
tree <- geiger::sim.bdtree(n = 10)
comm <- matrix(data = rpois(200, lambda = 1), 
               nrow = 20,
               ncol = 10,
               dimnames = list(paste("comm", 1:20, sep = "_"), tree$tip.label))

comm <- ifelse(comm >= 1, 1, 0)

# objetos -----------------------------------------------------------------

ada.obj <- 
  ada_core(x = comm, 
           phy = tree, 
           type = "continuous",
           lik.threshold = FALSE, 
           threshold = 0.6, 
           compute.node.by.sites = TRUE)

# calculo ada

calc_PD_ada(ada.obj = ada.obj, threshold = 0.6, comm.name = "comm_3")

PD_ada(ada.obj = ada.obj, thresold = 0.6, comm = comm)

