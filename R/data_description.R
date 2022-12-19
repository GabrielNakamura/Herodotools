#' Hypothetical phylogenetic tree 
#'
#' This is a made up data that simulates a phylogenetic tree in newick format. This tree contains
#' five species and the annotations in each node corresponds to the occurrence area of ancestors 
#'
#' @format A newick phylogenetic tree with 5 species
#
"toy_treeEx"

#' Phylogenetic tree in of Akodon genus
#' 
#' This is a phylogenetic tree containing 30 species from genus Akodon
#' 
#' @format newick phylogenetic tree
"akodon_newick"

#' Occurrence data for Akodon species
#' 
#' This is a occurrence data for species from Akodon genus
#' 
#' @format data frame containing 732 rows and 40 columns 
'akodon_sites'

#' Output from `evoregion` function with classification of 732 assemblages
#' 
#' Object containing a list with two objects:
#'     \describe{
#'     \item{$PCPS}{A list containing eigenvectors from PCPS analysis}
#'     \item{$Cluster_Evoregions}{A vector containing the groups of each assemblage}
#'     }
'regions'

#' Output from BioGeoBEARS
#' 
#' Object from a ancestral reconstruction model with BioGeoBEARS
#' 
#' @format A list obtained from BioGeoBEARS with Akodon species
'resDEC'

#' Phylogenetic tree with 285 species
#' 
#' Dated (my) phylogenetic tree in newick format containing 285 of sigmodontine species
#' 
#' @format Newick phylogenetic tree
'rodent_phylo'

#' Trait on foraging strata for small rodent
#' 
#' Character vector indicating the foraging strata for 214 species of small rodents
#' 
#' @format character vector with three categories:
#'     \describe{
#'         \item{Ar}{Arboreal}
#'         \item{G}{Ground level}
#'         \item{S}{Scansorial}
#'      }   
'trait'