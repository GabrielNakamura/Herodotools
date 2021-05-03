

#### R Script to implement Evoregions ####
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
### Define matrix P and Evoregions ###

# Calculate matrix P (containing phylogenetic composition for each cell/site) and PCPS
pcps.muro.bray <-
  PCPS::pcps(comm, phylodist = cophenetic(phy), method = "bray")
P <- pcps.muro.bray$P
values.bray <- pcps.muro.bray$values
# Define a threshold value (eigenvectors containing more than 5% of variation)
thresh.bray <- max(which(values.bray[, 2] >= 0.05))
cum.sum.thresh.bray <-
  cumsum(as.data.frame(values.bray[, 2]))[1:thresh.bray,][3]
vec.bray <- pcps.muro.bray$vectors

# First, find the maxinum number of clusters to be used in the analysis
# The max.n.clust chosen must have a high congruence (i.e. return always the same clusters between runs)
# Find max nclust
source("find.max.nclust.R")
find.max.number.cluster <-
  find.max.nclust(
    x = vec.bray[, 1:thresh.bray],
    threshold=thresh.bray,
    nperm = 1000,
    max.nclust = c(10, 15, 20, 25, 30),
    subset = 350,
    confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
  )

# Finding evoregions
clust.vec.bray <-
  adegenet::find.clusters(
    vec.bray[, 1:thresh.bray],
    clust = NULL,
    choose.n.clust = FALSE,
    n.pca = thresh.bray,
    method = "kmeans",
    stat = "BIC",
    n.iter = 1e7,
    criterion = "diffNgroup",
    max.n.clust = 10
  )
rownames(as.data.frame(clust.vec.bray$grp)) == rownames(esp)

#### FIGURE 1 ####
quartz()
plot(esp,
     col = clust.vec.bray$grp,
     pch = 15,
     cex = 1.2)
#### FIGURE 1 ####


###############################################################
# Define cell affiliation to evoregion (following Olivero et al. 2013)
groups.vec.bray <- clust.vec.bray$grp
G1 <- which(groups.vec.bray == 1)
G2 <- which(groups.vec.bray == 2)
G3 <- which(groups.vec.bray == 3)
G4 <- which(groups.vec.bray == 4)

dist.P.fuzzy <-
  as.matrix(vegan::vegdist(
    x = vec.bray,
    method = "euclidean",
    diag = T,
    upper = T
  ))

PG1 <- dist.P.fuzzy[G1, G1]
fuzzy.OGU.1 <-
  matrix(NA, nrow(PG1), 2, dimnames = list(rownames(PG1), c("fuzzy.belonging", "group")))
for (z in 1:nrow(PG1)) {
  dis <- as.data.frame(PG1[z,])[-z,]
  fuzzy.OGU.1[z, 1] <- mean(dis)
}

fuzzy.OGU.1.pad <- scales::rescale(fuzzy.OGU.1, c(0, 1))
fuzzy.OGU.1.pad[, 2] <- 1
hist(fuzzy.OGU.1.pad[, 1])

PG2 <- dist.P.fuzzy[G2, G2]
dim(PG2)
fuzzy.OGU.2 <-
  matrix(NA, nrow(PG2), 2, dimnames = list(rownames(PG2), c("fuzzy.belonging", "group")))
for (y in 1:nrow(PG2)) {
  dis <- as.data.frame(PG2[y,])[-y,]
  fuzzy.OGU.2[y, 1] <- mean(dis)
}

fuzzy.OGU.2.pad <- scales::rescale(fuzzy.OGU.2, c(0, 1))
fuzzy.OGU.2.pad[, 2] <- 2
hist(fuzzy.OGU.2.pad[, 1])

PG3 <- dist.P.fuzzy[G3, G3]
fuzzy.OGU.3 <-
  matrix(NA, nrow(PG3), 2, dimnames = list(rownames(PG3), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG3)) {
  dis <- as.data.frame(PG3[w,])[-w,]
  fuzzy.OGU.3[w, 1] <- mean(dis)
}

fuzzy.OGU.3.pad <- scales::rescale(fuzzy.OGU.3, c(0, 1))
fuzzy.OGU.3.pad[, 2] <- 3
hist(fuzzy.OGU.3.pad[, 1])

PG4 <- dist.P.fuzzy[G4, G4]
fuzzy.OGU.4 <-
  matrix(NA, nrow(PG4), 2, dimnames = list(rownames(PG4), c("fuzzy.belonging", "group")))
for (w in 1:nrow(PG4)) {
  dis <- as.data.frame(PG4[w,])[-w,]
  fuzzy.OGU.4[w, 1] <- mean(dis)
}

fuzzy.OGU.4.pad <- scales::rescale(fuzzy.OGU.4, c(0, 1))
fuzzy.OGU.4.pad[, 2] <- 4
hist(fuzzy.OGU.4.pad[, 1])

fuzzy.OGU.all.pad <-
  rbind(fuzzy.OGU.1.pad,
        fuzzy.OGU.2.pad,
        fuzzy.OGU.3.pad,
        fuzzy.OGU.4.pad)
org <-
  SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad, envir = esp)
fuzzy.OGU.all.pad.org <- org$community
fuzzy.OGU.all.pad.org[, 1] <- 1 - fuzzy.OGU.all.pad.org[, 1]
esp.fuzzy <- org$environmental
rownames(fuzzy.OGU.all.pad.org) == rownames(esp.fuzzy)

#### FIGURE 2 ####
quartz()
datCol <-
  topo.colors(10)[as.numeric(cut(fuzzy.OGU.all.pad.org[, 1], breaks = 10))]
plot(esp.fuzzy,
     col = datCol,
     pch = 15,
     cex = 1.0)
#### FIGURE 2 ####

summary(fuzzy.OGU.all.pad.org[, 1])

###############################################################
# Calculate association between species and evoregions
groups.vec.bray <- as.data.frame(groups.vec.bray)
n.groups <- length(clust.vec.bray$size)
dummy.groups.vec.bray <-
  matrix(NA, nrow = nrow(groups.vec.bray), ncol = n.groups)
rownames(dummy.groups.vec.bray) <- rownames(groups.vec.bray)
colnames(dummy.groups.vec.bray) <- 1:n.groups

for (i in 1:n.groups) {
  dummy.groups.vec.bray[, i] <- as.numeric(groups.vec.bray == i)
}

c <- t(comm) %*% dummy.groups.vec.bray
c.pad <- vegan::decostand(c, "total")
spp.groups.bray_0.6 <- ifelse(c.pad >= 0.6, 1, 0)
spp.groups.bray_0.7 <- ifelse(c.pad >= 0.7, 1, 0)
spp.groups.bray_0.8 <- ifelse(c.pad >= 0.8, 1, 0)
spp.groups.bray_0.9 <- ifelse(c.pad >= 0.9, 1, 0)

d_0.6 <- matrix(0, nrow(spp.groups.bray_0.6), ncol = 1)
rownames(d_0.6) <- rownames(spp.groups.bray_0.6)
for (k in 1:ncol(spp.groups.bray_0.6)) {
  d_0.6 <- ifelse(spp.groups.bray_0.6[, k] == 1, k, d_0.6)
}
d_0.6 <- ifelse(d_0.6 == 0, n.groups + 1, d_0.6)

d_0.7 <- matrix(0, nrow(spp.groups.bray_0.7), ncol = 1)
rownames(d_0.7) <- rownames(spp.groups.bray_0.7)
for (k in 1:ncol(spp.groups.bray_0.7)) {
  d_0.7 <- ifelse(spp.groups.bray_0.7[, k] == 1, k, d_0.7)
}
d_0.7 <- ifelse(d_0.7 == 0, n.groups + 1, d_0.7)

d_0.8 <- matrix(0, nrow(spp.groups.bray_0.8), ncol = 1)
rownames(d_0.8) <- rownames(spp.groups.bray_0.8)
for (k in 1:ncol(spp.groups.bray_0.8)) {
  d_0.8 <- ifelse(spp.groups.bray_0.8[, k] == 1, k, d_0.8)
}
d_0.8 <- ifelse(d_0.8 == 0, n.groups + 1, d_0.8)

d_0.9 <- matrix(0, nrow(spp.groups.bray_0.9), ncol = 1)
rownames(d_0.9) <- rownames(spp.groups.bray_0.9)
for (k in 1:ncol(spp.groups.bray_0.9)) {
  d_0.9 <- ifelse(spp.groups.bray_0.9[, k] == 1, k, d_0.9)
}
d_0.9 <- ifelse(d_0.9 == 0, n.groups + 1, d_0.9)

###############################################################
# Reconstructing ancestral affiliation using different thresholds
anc.tree_0.6 <-
  phytools::make.simmap(phy, d_0.6, model = "ER", nsim = 1000)
pd_0.6 <- summary(anc.tree_0.6)
plot(
  pd_0.6,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.6
)

anc.tree_0.7 <-
  phytools::make.simmap(phy, d_0.7, model = "ER", nsim = 1000)
pd_0.7 <- summary(anc.tree_0.7)

#### FIGURE 3 ####
quartz()
plot(
  pd_0.7,
  fsize = 0.2,
  cex = 0.2,
  ftype = "i",
  labels = d_0.7,
  type = "fan"
)
#### FIGURE 3 ####

anc.tree_0.8 <-
  phytools::make.simmap(phy, d_0.8, model = "ER", nsim = 1000)
pd_0.8 <- summary(anc.tree_0.8)
plot(
  pd_0.8,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.8
)

anc.tree_0.9 <-
  phytools::make.simmap(phy, d_0.9, model = "ER", nsim = 1000)
pd_0.9 <- summary(anc.tree_0.9)
plot(
  pd_0.9,
  fsize = 0.5,
  cex = 0.3,
  ftype = "i",
  labels = d_0.9
)