% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_spp_nodes.R
\name{spp_nodes}
\alias{spp_nodes}
\title{Species and their respective ancestral nodes}
\usage{
spp_nodes(tree, node.prefix = "N")
}
\arguments{
\item{tree}{Phylogenetic tree}

\item{node.prefix}{Character indicating the prefix to be used to name nodes, default is the letter "N"}
}
\value{
A matrix with species in lines and nodes in columns. 1 indicates that the node corresponds to the ancestor of that species and zero indicates that species
does not share the ancestor node
}
\description{
This function computes a matrix containing species and their respective ancestral nodes
}
