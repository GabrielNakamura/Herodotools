# Evoregions
Evoregions: mapping shifts in phylogenetic turnover across biogeographic regions. 
Published in [Maestri & Duarte 2020] (https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13492)

# Transition rates
The function implement three different tip-based metrics of trait evolution published in [Luza et al. 2021] (https://onlinelibrary.wiley.com/doi/10.1002/ece3.8476). The function allows to estimate a tip(species)-wise rate and time of trait transition across the branches of a phylogeny. The function currently makes use of the "make.simmap" algorithm ("phytools" package) to estimate trait change across the phylogeny branches. The number of changes, and the time in which these changes took place, are then used to calculate the three tip-based metrics, called "Transition Rates", "Stasis Time", and "Last Transition Time".


## Install

To install Evoregions package type:

`devtools::install_github("GabrielNakamura/Evoregions", ref = "main")`


