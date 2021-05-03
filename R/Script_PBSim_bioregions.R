

#### R Script to implement PBSim/PBSor and bioregions ####
#### for comparison with evoregions                   ####
#### An example using the world-wide distributed rats and mice (Muroidea) ####

###############################################################
# Load packages
library(ape)
library(picante)
library(PCPS)
library(adegenet)
library(vegan)
library(scales)
library(SYNCSA)
library(phytools)

###############################################################
### Load (i) geographic coordinates, (ii) spp x sites matrix, and (iii) phylogeny ###

## Load (i) geographic coordinates for each cell/site in the grid
# This is a two-column matrix with Longitude and Latitude for each cell
# Sites as rows and Long/Lat as columns
coord <- as.matrix(read.table("coords.R", h = T))

## Load (ii) spp x sites matrix
# This is a matrix of the species composition for each cell/site
# Sites as rows and species as columns
muro <- as.matrix(read.table("comm.R", h = T))
# Remove sites with fewer than three species to avoid artifacts
zero.comm <-
  which(rowSums(muro) <= 2) # remove cells with zero, one and two spp
muro.clean <- muro[-zero.comm, ]
dim(muro.clean)
esp <- coord[-zero.comm,]

## Load (iii) a phylogeny for the study group
tree.muro <- ape::read.tree("tree_muroidea_770sp.tre")
ape::is.ultrametric(tree.muro)
# Compare species in phylogeny with assemblage data and prune datasets to match one another
match <-
  picante::match.phylo.comm(tree.muro, muro.clean) #standardize species in phylo and comm
phy <- match$phy
comm <- match$comm
rownames(esp) == rownames(comm)
# After removing species because of the phylogeny, check again for cells with 0, 1, and 2 presences
zero.comm <-
  which(rowSums(comm) <= 2) # remove cells with fewer than three species
comm <- comm[-zero.comm, ]
esp <- esp[-zero.comm,]


###############################################################
# Calculate PBSim and defining phyloregions
source("beta.pd.decompo.R")
PBSim <-
  beta.pd.decompo(comm,
                  phy,
                  type = "PhyloSor",
                  output.dist = T,
                  random = F)
PBSim.mat <- PBSim$betadiv$PhyloSor_turn
PCPS.PBSim.wcmd <-
  PCPS::wcmdscale.org(
    PBSim.mat,
    squareroot = F,
    eig = TRUE,
    correlations = F
  )
values.PBSim <- PCPS.PBSim.wcmd$values
# Define a threshold value (eigenvectors containing more than 5% of variation)
thresh.PBSim <- max(which(values.PBSim[, 2] >= 0.05))
cum.sum.thresh.PBSim <-
  cumsum(as.data.frame(values.PBSim[, 2]))[1:thresh.PBSim,][3]
vec.PBSim <- PCPS.PBSim.wcmd$vectors
# Find max nclust
source("find.max.nclust.R")
find.max.number.cluster.PBSim <-
  find.max.nclust(
    x = vec.PBSim[, 1:thresh.PBSim],
    threshold=thresh.PBSim,
    nperm = 1000,
    max.nclust = c(10, 15, 20, 25, 30),
    subset = 350,
    confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
  )

# Finding phyloregions
clust.vec.PBSim <-
  adegenet::find.clusters(
    vec.PBSim[, 1:thresh.PBSim],
    clust = NULL,
    choose.n.clust = FALSE,
    n.pca = thresh.PBSim,
    method = "kmeans",
    stat = "BIC",
    n.iter = 1e7,
    criterion = "diffNgroup",
    max.n.clust = 10
  )
rownames(as.data.frame(clust.vec.PBSim$grp)) == rownames(esp)

#### FIGURE phyloregions PBSim ####
quartz()
plot(esp,
     col = clust.vec.PBSim$grp,
     pch = 15,
     cex = 1.2)
#### FIGURE phyloregions PBSim ####

###############################################################
# Define cell affiliation to phyloregion (following Olivero et al. 2013)
groups.vec.PBSim <- clust.vec.PBSim$grp
G1.PBSim <- which(groups.vec.PBSim == 1)
G2.PBSim <- which(groups.vec.PBSim == 2)
G3.PBSim <- which(groups.vec.PBSim == 3)
G4.PBSim <- which(groups.vec.PBSim == 4)
G5.PBSim <- which(groups.vec.PBSim == 5)
G6.PBSim <- which(groups.vec.PBSim == 6)


dist.P.fuzzy.PBSim <-
  as.matrix(vegan::vegdist(
    x = vec.PBSim,
    method = "euclidean",
    diag = T,
    upper = T
  ))

PG1.PBSim <- dist.P.fuzzy.PBSim[G1.PBSim, G1.PBSim]
fuzzy.OGU.1.PBSim <-
  matrix(NA, nrow(PG1.PBSim), 2, dimnames = list(rownames(PG1.PBSim), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG1.PBSim)) {
  dis <- as.data.frame(PG1.PBSim[z,])[-z,]
  fuzzy.OGU.1.PBSim[z, 1] <- mean(dis)
}

fuzzy.OGU.1.pad.PBSim <- scales::rescale(fuzzy.OGU.1.PBSim, c(0, 1))
fuzzy.OGU.1.pad.PBSim[, 2] <- 1
hist(fuzzy.OGU.1.pad.PBSim[, 1])

PG2.PBSim <- dist.P.fuzzy.PBSim[G2.PBSim, G2.PBSim]
dim(PG2.PBSim)
fuzzy.OGU.2.PBSim <-
  matrix(NA, nrow(PG2.PBSim), 2, dimnames = list(rownames(PG2.PBSim), c("fuzzy.belonging", "group")))
for (y in 1:nrow(PG2.PBSim)) {
  dis <- as.data.frame(PG2.PBSim[y,])[-y,]
  fuzzy.OGU.2.PBSim[y, 1] <- mean(dis)
}

fuzzy.OGU.2.pad.PBSim <- scales::rescale(fuzzy.OGU.2.PBSim, c(0, 1))
fuzzy.OGU.2.pad.PBSim[, 2] <- 2
hist(fuzzy.OGU.2.pad.PBSim[, 1])

PG3.PBSim <- dist.P.fuzzy.PBSim[G3.PBSim, G3.PBSim]
fuzzy.OGU.3.PBSim <-
  matrix(NA, nrow(PG3.PBSim), 2, dimnames = list(rownames(PG3.PBSim), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG3.PBSim)) {
  dis <- as.data.frame(PG3.PBSim[w,])[-w,]
  fuzzy.OGU.3.PBSim[w, 1] <- mean(dis)
}

fuzzy.OGU.3.pad.PBSim <- scales::rescale(fuzzy.OGU.3.PBSim, c(0, 1))
fuzzy.OGU.3.pad.PBSim[, 2] <- 3
hist(fuzzy.OGU.3.pad.PBSim[, 1])

PG4.PBSim <- dist.P.fuzzy.PBSim[G4.PBSim, G4.PBSim]
fuzzy.OGU.4.PBSim <-
  matrix(NA, nrow(PG4.PBSim), 2, dimnames = list(rownames(PG4.PBSim), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG4.PBSim)) {
  dis <- as.data.frame(PG4.PBSim[w,])[-w,]
  fuzzy.OGU.4.PBSim[w, 1] <- mean(dis)
}

fuzzy.OGU.4.pad.PBSim <- scales::rescale(fuzzy.OGU.4.PBSim, c(0, 1))
fuzzy.OGU.4.pad.PBSim[, 2] <- 4
hist(fuzzy.OGU.4.pad.PBSim[, 1])

PG5.PBSim <- dist.P.fuzzy.PBSim[G5.PBSim, G5.PBSim]
fuzzy.OGU.5.PBSim <-
  matrix(NA, nrow(PG5.PBSim), 2, dimnames = list(rownames(PG5.PBSim), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG5.PBSim)) {
  dis <- as.data.frame(PG5.PBSim[z,])[-z,]
  fuzzy.OGU.5.PBSim[z, 1] <- mean(dis)
}

fuzzy.OGU.5.pad.PBSim <- scales::rescale(fuzzy.OGU.5.PBSim, c(0, 1))
fuzzy.OGU.5.pad.PBSim[, 2] <- 5
hist(fuzzy.OGU.5.pad.PBSim[, 1])

PG6.PBSim <- dist.P.fuzzy.PBSim[G6.PBSim, G6.PBSim]
fuzzy.OGU.6.PBSim <-
  matrix(NA, nrow(PG6.PBSim), 2, dimnames = list(rownames(PG6.PBSim), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG6.PBSim)) {
  dis <- as.data.frame(PG6.PBSim[z,])[-z,]
  fuzzy.OGU.6.PBSim[z, 1] <- mean(dis)
}

fuzzy.OGU.6.pad.PBSim <- scales::rescale(fuzzy.OGU.6.PBSim, c(0, 1))
fuzzy.OGU.6.pad.PBSim[, 2] <- 6
hist(fuzzy.OGU.6.pad.PBSim[, 1])

fuzzy.OGU.all.pad.PBSim <-
  rbind(
    fuzzy.OGU.1.pad.PBSim,
    fuzzy.OGU.2.pad.PBSim,
    fuzzy.OGU.3.pad.PBSim,
    fuzzy.OGU.4.pad.PBSim,
    fuzzy.OGU.5.pad.PBSim,
    fuzzy.OGU.6.pad.PBSim
  )
org.PBSim <-
  SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad.PBSim, envir = esp)
fuzzy.OGU.all.pad.org.PBSim <- org.PBSim$community
fuzzy.OGU.all.pad.org.PBSim[, 1] <-
  1 - fuzzy.OGU.all.pad.org.PBSim[, 1]
esp.fuzzy.PBSim <- org.PBSim$environmental
rownames(fuzzy.OGU.all.pad.org.PBSim) == rownames(esp.fuzzy.PBSim)

#### FIGURE affiliation PBSim ####
quartz()
datCol <-
  topo.colors(10)[as.numeric(cut(fuzzy.OGU.all.pad.org.PBSim[, 1], breaks = 10))]
plot(esp.fuzzy.PBSim,
     col = datCol,
     pch = 15,
     cex = 1.0)
#### FIGURE affiliation PBSim ####

summary(fuzzy.OGU.all.pad.org.PBSim[, 1])

###############################################################
# Calculate association between species and phyloregions
groups.vec.PBSim <- as.data.frame(groups.vec.PBSim)
n.groups.PBSim <- length(clust.vec.PBSim$size)
dummy.groups.vec.PBSim <-
  matrix(NA, nrow = nrow(groups.vec.PBSim), ncol = n.groups.PBSim)
rownames(dummy.groups.vec.PBSim) <- rownames(groups.vec.PBSim)
colnames(dummy.groups.vec.PBSim) <- 1:n.groups.PBSim

for (i in 1:n.groups.PBSim) {
  dummy.groups.vec.PBSim[, i] <- as.numeric(groups.vec.PBSim == i)
}

c.PBSim <- t(comm) %*% dummy.groups.vec.PBSim
c.pad.PBSim <- vegan::decostand(c.PBSim, "total")
spp.groups.PBSim_0.6 <- ifelse(c.pad.PBSim >= 0.6, 1, 0)
spp.groups.PBSim_0.7 <- ifelse(c.pad.PBSim >= 0.7, 1, 0)
spp.groups.PBSim_0.8 <- ifelse(c.pad.PBSim >= 0.8, 1, 0)
spp.groups.PBSim_0.9 <- ifelse(c.pad.PBSim >= 0.9, 1, 0)

d_0.6.PBSim <- matrix(0, nrow(spp.groups.PBSim_0.6), ncol = 1)
rownames(d_0.6.PBSim) <- rownames(spp.groups.PBSim_0.6)
for (k in 1:ncol(spp.groups.PBSim_0.6)) {
  d_0.6.PBSim <-
    ifelse(spp.groups.PBSim_0.6[, k] == 1, k, d_0.6.PBSim)
}
d_0.6.PBSim <-
  ifelse(d_0.6.PBSim == 0, n.groups.PBSim + 1, d_0.6.PBSim)

d_0.7.PBSim <- matrix(0, nrow(spp.groups.PBSim_0.7), ncol = 1)
rownames(d_0.7.PBSim) <- rownames(spp.groups.PBSim_0.7)
for (k in 1:ncol(spp.groups.PBSim_0.7)) {
  d_0.7.PBSim <-
    ifelse(spp.groups.PBSim_0.7[, k] == 1, k, d_0.7.PBSim)
}
d_0.7.PBSim <-
  ifelse(d_0.7.PBSim == 0, n.groups.PBSim + 1, d_0.7.PBSim)

d_0.8.PBSim <- matrix(0, nrow(spp.groups.PBSim_0.8), ncol = 1)
rownames(d_0.8.PBSim) <- rownames(spp.groups.PBSim_0.8)
for (k in 1:ncol(spp.groups.PBSim_0.8)) {
  d_0.8.PBSim <- ifelse(spp.groups.PBSim_0.8[, k] == 1, k, d_0.8.PBSim)
}
d_0.8.PBSim <-
  ifelse(d_0.8.PBSim == 0, n.groups.PBSim + 1, d_0.8.PBSim)

d_0.9.PBSim <- matrix(0, nrow(spp.groups.PBSim_0.9), ncol = 1)
rownames(d_0.9.PBSim) <- rownames(spp.groups.PBSim_0.9)
for (k in 1:ncol(spp.groups.PBSim_0.9)) {
  d_0.9.PBSim <-
    ifelse(spp.groups.PBSim_0.9[, k] == 1, k, d_0.9.PBSim)
}
d_0.9.PBSim <-
  ifelse(d_0.9.PBSim == 0, n.groups.PBSim + 1, d_0.9.PBSim)

###############################################################
# Reconstructing ancestral affiliation using different thresholds
anc.tree_0.6.PBSim <-
  phytools::make.simmap(phy, d_0.6.PBSim, model = "ER", nsim = 1000)
pd_0.6.PBSim <- summary(anc.tree_0.6.PBSim)
plot(
  pd_0.6.PBSim,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.6.PBSim
)

anc.tree_0.7.PBSim <-
  phytools::make.simmap(phy, d_0.7.PBSim, model = "ER", nsim = 1000)
pd_0.7.PBSim <- summary(anc.tree_0.7.PBSim)

#### FIGURE Tree PBSim ####
quartz()
plot(
  pd_0.7.PBSim,
  fsize = 0.2,
  cex = 0.2,
  ftype = "i",
  labels = d_0.7.PBSim,
  type = "fan"
)
#### FIGURE Tree PBSim ####

anc.tree_0.8.PBSim <-
  phytools::make.simmap(phy, d_0.8.PBSim, model = "ER", nsim = 1000)
pd_0.8.PBSim <- summary(anc.tree_0.8.PBSim)
plot(
  pd_0.8.PBSim,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.8.PBSim
)

anc.tree_0.9.PBSim <-
  phytools::make.simmap(phy, d_0.9.PBSim, model = "ER", nsim = 1000)
pd_0.9.PBSim <- summary(anc.tree_0.9.PBSim)
plot(
  pd_0.9.PBSim,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.9.PBSim
)


###############################################################
# Calculate PBSor and defining phyloregions
PBSor.mat <- PBSim$betadiv$PhyloSor
PCPS.PBSor.wcmd <-
  PCPS::wcmdscale.org(
    PBSor.mat,
    squareroot = F,
    eig = TRUE,
    correlations = F
  )
values.PBSor <- PCPS.PBSor.wcmd$values
# Define a threshold value (eigenvectors containing more than 5% of variation)
thresh.PBSor <- max(which(values.PBSor[, 2] >= 0.05))
cum.sum.thresh.PBSor <-
  cumsum(as.data.frame(values.PBSor[, 2]))[1:thresh.PBSor,][3]
vec.PBSor <- PCPS.PBSor.wcmd$vectors
# Find max nclust
source("find.max.nclust.R")
find.max.number.cluster.PBSor <-
  find.max.nclust(
    x = vec.PBSor[, 1:thresh.PBSor],
    threshold=thresh.PBSor,
    nperm = 1000,
    max.nclust = c(10, 15, 20, 25, 30),
    subset = 350,
    confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
  )

# Finding phyloregions
clust.vec.PBSor <-
  adegenet::find.clusters(
    vec.PBSor[, 1:thresh.PBSor],
    clust = NULL,
    choose.n.clust = FALSE,
    n.pca = thresh.PBSor,
    method = "kmeans",
    stat = "BIC",
    n.iter = 1e7,
    criterion = "diffNgroup",
    max.n.clust = 10
  )
rownames(as.data.frame(clust.vec.PBSor$grp)) == rownames(esp)

#### FIGURE phyloregions PBSor ####
quartz()
plot(esp,
     col = clust.vec.PBSor$grp,
     pch = 15,
     cex = 1.2)
#### FIGURE phyloregions PBSor ####

###############################################################
# Define cell affiliation to phyloregion (following Olivero et al. 2013)
groups.vec.PBSor <- clust.vec.PBSor$grp
G1.PBSor <- which(groups.vec.PBSor == 1)
G2.PBSor <- which(groups.vec.PBSor == 2)
G3.PBSor <- which(groups.vec.PBSor == 3)
G4.PBSor <- which(groups.vec.PBSor == 4)
G5.PBSor <- which(groups.vec.PBSor == 5)
G6.PBSor <- which(groups.vec.PBSor == 6)
G7.PBSor <- which(groups.vec.PBSor == 7)

dist.P.fuzzy.PBSor <-
  as.matrix(vegan::vegdist(
    x = vec.PBSor,
    method = "euclidean",
    diag = T,
    upper = T
  ))

PG1.PBSor <- dist.P.fuzzy.PBSor[G1.PBSor, G1.PBSor]
fuzzy.OGU.1.PBSor <-
  matrix(NA, nrow(PG1.PBSor), 2, dimnames = list(rownames(PG1.PBSor), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG1.PBSor)) {
  dis <- as.data.frame(PG1.PBSor[z,])[-z,]
  fuzzy.OGU.1.PBSor[z, 1] <- mean(dis)
}

fuzzy.OGU.1.pad.PBSor <- scales::rescale(fuzzy.OGU.1.PBSor, c(0, 1))
fuzzy.OGU.1.pad.PBSor[, 2] <- 1
hist(fuzzy.OGU.1.pad.PBSor[, 1])

PG2.PBSor <- dist.P.fuzzy.PBSor[G2.PBSor, G2.PBSor]
dim(PG2.PBSor)
fuzzy.OGU.2.PBSor <-
  matrix(NA, nrow(PG2.PBSor), 2, dimnames = list(rownames(PG2.PBSor), c("fuzzy.belonging", "group")))
for (y in 1:nrow(PG2.PBSor)) {
  dis <- as.data.frame(PG2.PBSor[y,])[-y,]
  fuzzy.OGU.2.PBSor[y, 1] <- mean(dis)
}

fuzzy.OGU.2.pad.PBSor <- scales::rescale(fuzzy.OGU.2.PBSor, c(0, 1))
fuzzy.OGU.2.pad.PBSor[, 2] <- 2
hist(fuzzy.OGU.2.pad.PBSor[, 1])

PG3.PBSor <- dist.P.fuzzy.PBSor[G3.PBSor, G3.PBSor]
fuzzy.OGU.3.PBSor <-
  matrix(NA, nrow(PG3.PBSor), 2, dimnames = list(rownames(PG3.PBSor), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG3.PBSor)) {
  dis <- as.data.frame(PG3.PBSor[w,])[-w,]
  fuzzy.OGU.3.PBSor[w, 1] <- mean(dis)
}

fuzzy.OGU.3.pad.PBSor <- scales::rescale(fuzzy.OGU.3.PBSor, c(0, 1))
fuzzy.OGU.3.pad.PBSor[, 2] <- 3
hist(fuzzy.OGU.3.pad.PBSor[, 1])

PG4.PBSor <- dist.P.fuzzy.PBSor[G4.PBSor, G4.PBSor]
fuzzy.OGU.4.PBSor <-
  matrix(NA, nrow(PG4.PBSor), 2, dimnames = list(rownames(PG4.PBSor), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG4.PBSor)) {
  dis <- as.data.frame(PG4.PBSor[w,])[-w,]
  fuzzy.OGU.4.PBSor[w, 1] <- mean(dis)
}

fuzzy.OGU.4.pad.PBSor <- scales::rescale(fuzzy.OGU.4.PBSor, c(0, 1))
fuzzy.OGU.4.pad.PBSor[, 2] <- 4
hist(fuzzy.OGU.4.pad.PBSor[, 1])

PG5.PBSor <- dist.P.fuzzy.PBSor[G5.PBSor, G5.PBSor]
fuzzy.OGU.5.PBSor <-
  matrix(NA, nrow(PG5.PBSor), 2, dimnames = list(rownames(PG5.PBSor), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG5.PBSor)) {
  dis <- as.data.frame(PG5.PBSor[z,])[-z,]
  fuzzy.OGU.5.PBSor[z, 1] <- mean(dis)
}

fuzzy.OGU.5.pad.PBSor <- scales::rescale(fuzzy.OGU.5.PBSor, c(0, 1))
fuzzy.OGU.5.pad.PBSor[, 2] <- 5
hist(fuzzy.OGU.5.pad.PBSor[, 1])

PG6.PBSor <- dist.P.fuzzy.PBSor[G6.PBSor, G6.PBSor]
fuzzy.OGU.6.PBSor <-
  matrix(NA, nrow(PG6.PBSor), 2, dimnames = list(rownames(PG6.PBSor), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG6.PBSor)) {
  dis <- as.data.frame(PG6.PBSor[z,])[-z,]
  fuzzy.OGU.6.PBSor[z, 1] <- mean(dis)
}

fuzzy.OGU.6.pad.PBSor <- scales::rescale(fuzzy.OGU.6.PBSor, c(0, 1))
fuzzy.OGU.6.pad.PBSor[, 2] <- 6
hist(fuzzy.OGU.6.pad.PBSor[, 1])

PG7.PBSor <- dist.P.fuzzy.PBSor[G7.PBSor, G7.PBSor]
fuzzy.OGU.7.PBSor <-
  matrix(NA, nrow(PG7.PBSor), 2, dimnames = list(rownames(PG7.PBSor), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG7.PBSor)) {
  dis <- as.data.frame(PG7.PBSor[z,])[-z,]
  fuzzy.OGU.7.PBSor[z, 1] <- mean(dis)
}

fuzzy.OGU.7.pad.PBSor <- scales::rescale(fuzzy.OGU.7.PBSor, c(0, 1))
fuzzy.OGU.7.pad.PBSor[, 2] <- 7
hist(fuzzy.OGU.7.pad.PBSor[, 1])

fuzzy.OGU.all.pad.PBSor <-
  rbind(
    fuzzy.OGU.1.pad.PBSor,
    fuzzy.OGU.2.pad.PBSor,
    fuzzy.OGU.3.pad.PBSor,
    fuzzy.OGU.4.pad.PBSor,
    fuzzy.OGU.5.pad.PBSor,
    fuzzy.OGU.6.pad.PBSor,
    fuzzy.OGU.7.pad.PBSor
  )
org.PBSor <-
  SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad.PBSor, envir = esp)
fuzzy.OGU.all.pad.org.PBSor <- org.PBSor$community
fuzzy.OGU.all.pad.org.PBSor[, 1] <-
  1 - fuzzy.OGU.all.pad.org.PBSor[, 1]
esp.fuzzy.PBSor <- org.PBSor$environmental
rownames(fuzzy.OGU.all.pad.org.PBSor) == rownames(esp.fuzzy.PBSor)

#### FIGURE affiliation PBSor ####
quartz()
datCol <-
  topo.colors(10)[as.numeric(cut(fuzzy.OGU.all.pad.org.PBSor[, 1], breaks = 10))]
plot(esp.fuzzy.PBSor,
     col = datCol,
     pch = 15,
     cex = 1.0)
#### FIGURE affiliation PBSor ####

summary(fuzzy.OGU.all.pad.org.PBSor[, 1])

###############################################################
# Calculate association between species and phyloregions
groups.vec.PBSor <- as.data.frame(groups.vec.PBSor)
n.groups.PBSor <- length(clust.vec.PBSor$size)
dummy.groups.vec.PBSor <-
  matrix(NA, nrow = nrow(groups.vec.PBSor), ncol = n.groups.PBSor)
rownames(dummy.groups.vec.PBSor) <- rownames(groups.vec.PBSor)
colnames(dummy.groups.vec.PBSor) <- 1:n.groups.PBSor

for (i in 1:n.groups.PBSor) {
  dummy.groups.vec.PBSor[, i] <- as.numeric(groups.vec.PBSor == i)
}

c.PBSor <- t(comm) %*% dummy.groups.vec.PBSor
c.pad.PBSor <- vegan::decostand(c.PBSor, "total")
spp.groups.PBSor_0.6 <- ifelse(c.pad.PBSor >= 0.6, 1, 0)
spp.groups.PBSor_0.7 <- ifelse(c.pad.PBSor >= 0.7, 1, 0)
spp.groups.PBSor_0.8 <- ifelse(c.pad.PBSor >= 0.8, 1, 0)
spp.groups.PBSor_0.9 <- ifelse(c.pad.PBSor >= 0.9, 1, 0)

d_0.6.PBSor <- matrix(0, nrow(spp.groups.PBSor_0.6), ncol = 1)
rownames(d_0.6.PBSor) <- rownames(spp.groups.PBSor_0.6)
for (k in 1:ncol(spp.groups.PBSor_0.6)) {
  d_0.6.PBSor <-
    ifelse(spp.groups.PBSor_0.6[, k] == 1, k, d_0.6.PBSor)
}
d_0.6.PBSor <-
  ifelse(d_0.6.PBSor == 0, n.groups.PBSor + 1, d_0.6.PBSor)

d_0.7.PBSor <- matrix(0, nrow(spp.groups.PBSor_0.7), ncol = 1)
rownames(d_0.7.PBSor) <- rownames(spp.groups.PBSor_0.7)
for (k in 1:ncol(spp.groups.PBSor_0.7)) {
  d_0.7.PBSor <-
    ifelse(spp.groups.PBSor_0.7[, k] == 1, k, d_0.7.PBSor)
}
d_0.7.PBSor <-
  ifelse(d_0.7.PBSor == 0, n.groups.PBSor + 1, d_0.7.PBSor)

d_0.8.PBSor <- matrix(0, nrow(spp.groups.PBSor_0.8), ncol = 1)
rownames(d_0.8.PBSor) <- rownames(spp.groups.PBSor_0.8)
for (k in 1:ncol(spp.groups.PBSor_0.8)) {
  d_0.8.PBSor <- ifelse(spp.groups.PBSor_0.8[, k] == 1, k, d_0.8.PBSor)
}
d_0.8.PBSor <-
  ifelse(d_0.8.PBSor == 0, n.groups.PBSor + 1, d_0.8.PBSor)

d_0.9.PBSor <- matrix(0, nrow(spp.groups.PBSor_0.9), ncol = 1)
rownames(d_0.9.PBSor) <- rownames(spp.groups.PBSor_0.9)
for (k in 1:ncol(spp.groups.PBSor_0.9)) {
  d_0.9.PBSor <-
    ifelse(spp.groups.PBSor_0.9[, k] == 1, k, d_0.9.PBSor)
}
d_0.9.PBSor <-
  ifelse(d_0.9.PBSor == 0, n.groups.PBSor + 1, d_0.9.PBSor)

###############################################################
# Reconstructing ancestral affiliation using different thresholds
anc.tree_0.6.PBSor <-
  phytools::make.simmap(phy, d_0.6.PBSor, model = "ER", nsim = 1000)
pd_0.6.PBSor <- summary(anc.tree_0.6.PBSor)
plot(
  pd_0.6.PBSor,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.6.PBSor
)

anc.tree_0.7.PBSor <-
  phytools::make.simmap(phy, d_0.7.PBSor, model = "ER", nsim = 1000)
pd_0.7.PBSor <- summary(anc.tree_0.7.PBSor)

#### FIGURE Tree PBSor ####
quartz()
plot(
  pd_0.7.PBSor,
  fsize = 0.2,
  cex = 0.2,
  ftype = "i",
  labels = d_0.7.PBSor,
  type = "fan"
)
#### FIGURE Tree PBSor ####

anc.tree_0.8.PBSor <-
  phytools::make.simmap(phy, d_0.8.PBSor, model = "ER", nsim = 1000)
pd_0.8.PBSor <- summary(anc.tree_0.8.PBSor)
plot(
  pd_0.8.PBSor,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.8.PBSor
)

anc.tree_0.9.PBSor <-
  phytools::make.simmap(phy, d_0.9.PBSor, model = "ER", nsim = 1000)
pd_0.9.PBSor <- summary(anc.tree_0.9.PBSor)
plot(
  pd_0.9.PBSor,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.9.PBSor
)


###############################################################
### Find bioregions using species composition alone BetaSOR ###
pcW.comm <-
  vegan::capscale(comm ~ 1, distance = "bray", sqrt.dist = T)
values.comm <- vegan::eigenvals(pcW.comm)
vectors.comm <-
  vegan::scores(pcW.comm, choices = c(1:length(values.comm)))$sites
rel.val.comm <- round(values.comm / sum(values.comm), digits = 5)
cum.val.comm <-
  round(cumsum(values.comm) / sum(values.comm), digits = 5)

thresh.pcw <- max(which(rel.val.comm >= 0.05))

# Find max nclust
source("find.max.nclust.R")
find.max.number.cluster.bioregion <-
    find.max.nclust(
    x = vectors.comm,
    threshold = thresh.pcW,
    nperm = 1000,
    max.nclust = c(10, 15, 20, 30, 50),
    subset = 350,
    confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
  )

clust.vec.comm <-
  adegenet::find.clusters(
    vectors.comm[, 1:thresh.pcw],
    clust = NULL,
    choose.n.clust = FALSE,
    n.pca = thresh.pcw,
    method = "kmeans",
    stat = "BIC",
    n.iter = 1e7,
    criterion = "diffNgroup",
    max.n.clust = 20
  )
rownames(as.data.frame(clust.vec.comm$grp)) == rownames(esp)

#### FIGURE Bioregions BetaSOR ####
datCol <-
  colorRampPalette(rainbow(length(clust.vec.comm$size), alpha = 0.1))(length(clust.vec.comm$size))[clust.vec.comm$grp]
plot(esp, col = datCol, pch = 15, cex = 1.2)
#### FIGURE Bioregions BetaSOR ####


###############################################################
# Define cell affiliation to bioregion
groups.vec.comm <- clust.vec.comm$grp
G1.comm <- which(groups.vec.comm == 1)
G2.comm <- which(groups.vec.comm == 2)
G3.comm <- which(groups.vec.comm == 3)
G4.comm <- which(groups.vec.comm == 4)

dist.P.comm <-
  as.matrix(vegan::vegdist(
    x = vectors.comm,
    method = "euclidean",
    diag = T,
    upper = T
  ))

PG1.comm <- dist.P.comm[G1.comm, G1.comm]
comm.OGU.1 <-
  matrix(NA, nrow(PG1.comm), 2, dimnames = list(rownames(PG1.comm), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG1.comm)) {
  dis <- as.data.frame(PG1.comm[z,])[-z,]
  comm.OGU.1[z, 1] <- mean(dis)
}

comm.OGU.1.pad <- scales::rescale(comm.OGU.1, c(0, 1))
comm.OGU.1.pad[, 2] <- 1
hist(comm.OGU.1.pad[, 1])

PG2.comm <- dist.P.comm[G2.comm, G2.comm]
dim(PG2.comm)
comm.OGU.2 <-
  matrix(NA, nrow(PG2.comm), 2, dimnames = list(rownames(PG2.comm), c("fuzzy.belonging", "group")))
for (y in 1:nrow(PG2.comm)) {
  dis <- as.data.frame(PG2.comm[y,])[-y,]
  comm.OGU.2[y, 1] <- mean(dis)
}

comm.OGU.2.pad <- scales::rescale(comm.OGU.2, c(0, 1))
comm.OGU.2.pad[, 2] <- 2
hist(comm.OGU.2.pad[, 1])

PG3.comm <- dist.P.comm[G3.comm, G3.comm]
comm.OGU.3 <-
  matrix(NA, nrow(PG3.comm), 2, dimnames = list(rownames(PG3.comm), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG3.comm)) {
  dis <- as.data.frame(PG3.comm[w,])[-w,]
  comm.OGU.3[w, 1] <- mean(dis)
}

comm.OGU.3.pad <- scales::rescale(comm.OGU.3, c(0, 1))
comm.OGU.3.pad[, 2] <- 3
hist(comm.OGU.3.pad[, 1])

PG4.comm <- dist.P.comm[G4.comm, G4.comm]
comm.OGU.4 <-
  matrix(NA, nrow(PG4.comm), 2, dimnames = list(rownames(PG4.comm), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG4.comm)) {
  dis <- as.data.frame(PG4.comm[w,])[-w,]
  comm.OGU.4[w, 1] <- mean(dis)
}

comm.OGU.4.pad <- scales::rescale(comm.OGU.4, c(0, 1))
comm.OGU.4.pad[, 2] <- 4
hist(comm.OGU.4.pad[, 1])

comm.OGU.all.pad <-
  rbind(comm.OGU.1.pad,
        comm.OGU.2.pad,
        comm.OGU.3.pad,
        comm.OGU.4.pad)
org <-
  SYNCSA::organize.syncsa(comm = comm.OGU.all.pad, envir = esp)
comm.OGU.all.pad.org <- org$community
comm.OGU.all.pad.org[, 1] <- 1 - comm.OGU.all.pad.org[, 1]
esp.comm <- org$environmental
rownames(comm.OGU.all.pad.org) == rownames(esp.comm)

#### FIGURE affiliation Bioregions BetaSOR ####
quartz()
datCol <-
  topo.colors(10)[as.numeric(cut(comm.OGU.all.pad.org[, 1], breaks = 10))]
plot(esp.comm,
     col = datCol,
     pch = 15,
     cex = 1.0)
#### FIGURE affiliation Bioregions BetaSOR ####

summary(comm.OGU.all.pad.org[, 1])


###############################################################
# Calculate association between species and bioregions
groups.vec.comm <- as.data.frame(clust.vec.comm$grp)
n.groups.comm <- length(clust.vec.comm$size)
dummy.groups.vec.comm <-
  matrix(NA, nrow = nrow(groups.vec.comm), ncol = n.groups.comm)
rownames(dummy.groups.vec.comm) <- rownames(groups.vec.comm)
colnames(dummy.groups.vec.comm) <- 1:n.groups.comm

for (i in 1:n.groups.comm) {
  dummy.groups.vec.comm[, i] <- as.numeric(groups.vec.comm == i)
}

c.comm <- t(comm) %*% dummy.groups.vec.comm
c.pad.comm <- vegan::decostand(c.comm, "total")
spp.groups.comm_0.7 <- ifelse(c.pad.comm >= 0.7, 1, 0) # 70%

d_0.7.comm <- matrix(0, nrow(spp.groups.comm_0.7), ncol = 1)
rownames(d_0.7.comm) <- rownames(spp.groups.comm_0.7)
for (k in 1:ncol(spp.groups.comm_0.7)) {
  d_0.7.comm <- ifelse(spp.groups.comm_0.7[, k] == 1, k, d_0.7.comm)
}
d_0.7.comm <- ifelse(d_0.7.comm == 0, n.groups + 1, d_0.7.comm)

anc.tree_0.7.comm <-
  phytools::make.simmap(phy, d_0.7.comm, model = "ER", nsim = 1000)
pd_0.7.comm <- summary(anc.tree_0.7.comm)

#### FIGURE Tree bioregions BetaSOR ####
quartz()
plot(
  pd_0.7.comm,
  fsize = 0.2,
  cex = 0.2,
  ftype = "i",
  labels = d_0.7.comm,
  type = "fan"
)
#### FIGURE Tree bioregions BetaSOR ####


###############################################################
### Find bioregions using species composition alone betaSIM ###
require(betapart)
Simpson<-beta.pair(comm,index.family="sorensen")
Sim<-Simpson$beta.sim

PCPS.Sim.wcmd<-PCPS::wcmdscale.org(Sim,squareroot=F,eig=TRUE,correlations=F)
values.Sim<-PCPS.Sim.wcmd$values
# Define a threshold value (eigenvectors containing more than 5% of variation)
thresh.Sim<- max(which(values.Sim[, 2] >= 0.05))
cum.sum.thresh.Sim <-cumsum(as.data.frame(values.Sim[, 2]))[1:thresh.Sim, ][3]
vec.Sim<-PCPS.Sim.wcmd$vectors

# Find max nclust
source("find.max.nclust.R")

find.max.number.cluster.Sim.350<-find.max.nclust(x=vec.Sim,threshold=thresh.Sim,
                                                 runs=1000,max.nclust=c(10,15,20,25,30),subset=350,
                                                 confidence.level=c(0.7,0.8,0.9,0.95,0.99))

clust.vec.Sim <-
  adegenet::find.clusters(
    vec.Sim[, 1:thresh.Sim],
    clust = NULL,
    choose.n.clust = FALSE,
    n.pca = thresh.Sim,
    method = "kmeans",
    stat = "BIC",
    n.iter = 1e7,
    criterion = "diffNgroup",
    max.n.clust = 10
  )
rownames(as.data.frame(clust.vec.Sim$grp)) == rownames(esp)

#### FIGURE Bioregions BetaSIM ####
datCol <-
  colorRampPalette(rainbow(length(clust.vec.Sim$size), alpha = 0.1))(length(clust.vec.Sim$size))[clust.vec.Sim$grp]
plot(esp, col = datCol, pch = 15, cex = 1.2)
#### FIGURE Bioregions BetaSIM ####

###############################################################
# Define cell affiliation to evoregion (following Olivero et al. 2013)
groups.vec.Sim <- clust.vec.Sim$grp
G1.Sim <- which(groups.vec.Sim == 1)
G2.Sim <- which(groups.vec.Sim == 2)
G3.Sim <- which(groups.vec.Sim == 3)
G4.Sim <- which(groups.vec.Sim == 4)
G5.Sim <- which(groups.vec.Sim == 5)
G6.Sim <- which(groups.vec.Sim == 6)

dist.Sim.fuzzy <-
  as.matrix(vegan::vegdist(
    x = vec.Sim,
    method = "euclidean",
    diag = T,
    upper = T
  ))

PG1.Sim <- dist.Sim.fuzzy[G1.Sim, G1.Sim]
fuzzy.OGU.1.Sim <-
  matrix(NA, nrow(PG1.Sim), 2, dimnames = list(rownames(PG1.Sim), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG1.Sim)) {
  dis <- as.data.frame(PG1.Sim[z, ])[-z, ]
  fuzzy.OGU.1.Sim[z, 1] <- mean(dis)
}

fuzzy.OGU.1.pad.Sim <- scales::rescale(fuzzy.OGU.1.Sim, c(0, 1))
fuzzy.OGU.1.pad.Sim[, 2] <- 1
hist(fuzzy.OGU.1.pad.Sim[, 1])

PG2.Sim <- dist.Sim.fuzzy[G2.Sim, G2.Sim]
dim(PG2.Sim)
fuzzy.OGU.2.Sim <-
  matrix(NA, nrow(PG2.Sim), 2, dimnames = list(rownames(PG2.Sim), c("fuzzy.belonging", "group")))
for (y in 1:nrow(PG2.Sim)) {
  dis <- as.data.frame(PG2.Sim[y, ])[-y, ]
  fuzzy.OGU.2.Sim[y, 1] <- mean(dis)
}

fuzzy.OGU.2.pad.Sim <- scales::rescale(fuzzy.OGU.2.Sim, c(0, 1))
fuzzy.OGU.2.pad.Sim[, 2] <- 2
hist(fuzzy.OGU.2.pad.Sim[, 1])

PG3.Sim <- dist.Sim.fuzzy[G3.Sim, G3.Sim]
fuzzy.OGU.3.Sim <-
  matrix(NA, nrow(PG3.Sim), 2, dimnames = list(rownames(PG3.Sim), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG3.Sim)) {
  dis <- as.data.frame(PG3.Sim[w, ])[-w, ]
  fuzzy.OGU.3.Sim[w, 1] <- mean(dis)
}

fuzzy.OGU.3.pad.Sim <- scales::rescale(fuzzy.OGU.3.Sim, c(0, 1))
fuzzy.OGU.3.pad.Sim[, 2] <- 3
hist(fuzzy.OGU.3.pad.Sim[, 1])

PG4.Sim <- dist.Sim.fuzzy[G4.Sim, G4.Sim]
fuzzy.OGU.4.Sim <-
  matrix(NA, nrow(PG4.Sim), 2, dimnames = list(rownames(PG4.Sim), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG4.Sim)) {
  dis <- as.data.frame(PG4.Sim[w, ])[-w, ]
  fuzzy.OGU.4.Sim[w, 1] <- mean(dis)
}

fuzzy.OGU.4.pad.Sim <- scales::rescale(fuzzy.OGU.4.Sim, c(0, 1))
fuzzy.OGU.4.pad.Sim[, 2] <- 4
hist(fuzzy.OGU.4.pad.Sim[, 1])

PG5.Sim <- dist.Sim.fuzzy[G5.Sim, G5.Sim]
fuzzy.OGU.5.Sim <-
  matrix(NA, nrow(PG5.Sim), 2, dimnames = list(rownames(PG5.Sim), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG5.Sim)) {
  dis <- as.data.frame(PG5.Sim[w, ])[-w, ]
  fuzzy.OGU.5.Sim[w, 1] <- mean(dis)
}

fuzzy.OGU.5.pad.Sim <- scales::rescale(fuzzy.OGU.5.Sim, c(0, 1))
fuzzy.OGU.5.pad.Sim[, 2] <- 5
hist(fuzzy.OGU.5.pad.Sim[, 1])

PG6.Sim <- dist.Sim.fuzzy[G6.Sim, G6.Sim]
fuzzy.OGU.6.Sim <-
  matrix(NA, nrow(PG6.Sim), 2, dimnames = list(rownames(PG6.Sim), c("fuzzy.belonging", "group")))
for (s in 1:nrow(PG6.Sim)) {
  dis <- as.data.frame(PG6.Sim[s, ])[-s, ]
  fuzzy.OGU.6.Sim[s, 1] <- mean(dis)
}

fuzzy.OGU.6.pad.Sim <- scales::rescale(fuzzy.OGU.6.Sim, c(0, 1))
fuzzy.OGU.6.pad.Sim[, 2] <- 6
hist(fuzzy.OGU.6.pad.Sim[, 1])


fuzzy.OGU.all.pad.Sim <-
  rbind(fuzzy.OGU.1.pad.Sim,
        fuzzy.OGU.2.pad.Sim,
        fuzzy.OGU.3.pad.Sim,
        fuzzy.OGU.4.pad.Sim,
        fuzzy.OGU.5.pad.Sim,
        fuzzy.OGU.6.pad.Sim
  )
org.Sim <- SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad.Sim, envir = esp)
fuzzy.OGU.all.pad.org.Sim <- org.Sim$community
fuzzy.OGU.all.pad.org.Sim[, 1] <- 1 - fuzzy.OGU.all.pad.org.Sim[, 1]
esp.fuzzy.Sim <- org.Sim$environmental
rownames(fuzzy.OGU.all.pad.org.Sim) == rownames(esp.fuzzy.Sim)

#### FIGURE affiliation Bioregions BetaSIM ####
quartz()
datCol <-
  topo.colors(10)[as.numeric(cut(fuzzy.OGU.all.pad.org.Sim[, 1], breaks = 10))]
plot(esp.fuzzy.Sim,
     col = datCol,
     pch = 15,
     cex = 1.0)
#### FIGURE affiliation Bioregions BetaSIM ####

###############################################################
# Calculate association between species and bioregions BetaSIM
groups.vec.comm <- as.data.frame(clust.vec.Sim$grp)
n.groups.comm <- length(clust.vec.Sim$size)
dummy.groups.vec.comm <-
  matrix(NA, nrow = nrow(groups.vec.comm), ncol = n.groups.comm)
rownames(dummy.groups.vec.comm) <- rownames(groups.vec.comm)
colnames(dummy.groups.vec.comm) <- 1:n.groups.comm

for (i in 1:n.groups.comm) {
  dummy.groups.vec.comm[, i] <- as.numeric(groups.vec.comm == i)
}

c.comm <- t(comm) %*% dummy.groups.vec.comm
c.pad.comm <- vegan::decostand(c.comm, "total")
spp.groups.comm_0.7 <- ifelse(c.pad.comm >= 0.7, 1, 0) # 70%

d_0.7.comm <- matrix(0, nrow(spp.groups.comm_0.7), ncol = 1)
rownames(d_0.7.comm) <- rownames(spp.groups.comm_0.7)
for (k in 1:ncol(spp.groups.comm_0.7)) {
  d_0.7.comm <- ifelse(spp.groups.comm_0.7[, k] == 1, k, d_0.7.comm)
}
d_0.7.comm <- ifelse(d_0.7.comm == 0, n.groups.comm + 1, d_0.7.comm)

anc.tree_0.7.comm <-
  phytools::make.simmap(phy, d_0.7.comm, model = "ER", nsim = 1000)
pd_0.7.comm <- summary(anc.tree_0.7.comm)

#### FIGURE Tree bioregions BetaSIM ####
quartz()
plot(
  pd_0.7.comm,
  fsize = 0.2,
  cex = 0.2,
  ftype = "i",
  labels = d_0.7.comm,
  type = "fan"
)
#### FIGURE Tree bioregions BetaSIM ####
