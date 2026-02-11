# Global variables used via NSE (dplyr / ggplot2 / base)
# Declared to satisfy R CMD check
# See: Writing R Extensions, "Non-standard evaluation"

utils::globalVariables(
  c(
    # tidy evaluation / data masking
    ".data",
    
    # base / table conversions
    "Freq",
    
    # spatial variables
    "LAT", "LONG", "grids",
    
    # phylogeny / tree structure
    "parent", "node", "child", "branch.length", "label",
    "nodenum_at_top_of_branch", "nodepath",
    
    # biogeography / dispersal
    "area", "pres", "dispersal", "new_rangetxt",
    
    # PD partitioning
    "ancestor", "descendant", "ancestor1", "descendent1",
    "partition.IS", "partition.IM", "partition.EM", "partition.ESD",
    
    # temporal / event data
    "event_time", "event_type", "mean_age_arrival", "unit"
 
  )
)
