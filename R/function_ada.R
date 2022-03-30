#' Ancestral Diversity Analysis
#'
#' @param x Community occurrence matrix. Rows are sites and columns are species
#' @param phy Phylogenetic tree
#' @param sp.bin Character indicating the methods to be used to compute the number of time slices in which metrics will be computed. Default is "Sturges"
#' @param marginal Logical indicating 
#' @param lik.threshold Logical indicating if some threshold value will be used to select the likelihood values of ancestor species in a site. Default is TRUE
#' @param threshold Scalar indicating the threshold value used to select the likelihood of a given species be present at a site
#' @param compute.fields Logical, indicates if phylogenetic fields will be computed. Default is FALSE
#'
#' @return 
#' @export
#'
#' @examples
#' 
ada <- function(x,
                phy, 
                sp.bin = "Sturges", # eu mudaria para breaks aqui 
                marginal = FALSE, 
                lik.threshold = FALSE, 
                threshold = 0.7, 
                compute.fields = F, 
                plot.results = TRUE){
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
  
  #h <- hist(ages,breaks=round(length(ages)/sp.bin,0),plot=F)
  h <- hist(ages, breaks = sp.bin, plot = F)
  #hist(ages, breaks = 20)
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
}
  
  