% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_ancestral_state.R
\name{calc_ancestral_state}
\alias{calc_ancestral_state}
\title{Ancestral State per Assemblage}
\usage{
calc_ancestral_state(tree, ancestral.area, prefix = "N")
}
\arguments{
\item{tree}{A newick phylogenetic tree object}

\item{ancestral.area}{One column data frame. Lines are nodes the column are the biomes/region of occurrence
for each ancestors. Can be obtained by using (\code{\link{get_node_range_BioGeoBEARS}})}

\item{prefix}{A single character string to be used to name nodes}
}
\value{
A data frame with assemblages in lines and nodes in columns.
Each cell contains the ancestral area/Ecoregion of occurrence for each node and its respective species.
}
\description{
Ancestral State per Assemblage
}
\examples{
biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
calc_ancestral_state(toy_treeEx, ancestral_area_toy)
}
\author{
Gabriel Nakamura \href{mailto:gabriel.nakamura.souza@gmail.com}{gabriel.nakamura.souza@gmail.com}
}
