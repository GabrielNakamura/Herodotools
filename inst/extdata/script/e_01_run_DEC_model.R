#######################################################
# SETUP -- libraries/BioGeoBEARS updates
#######################################################

# Load the package (after installation, see above).
library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(dplyr)
library(here)


#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################

phy.path <-  here("inst", "extdata", "akodon.new")

trfn = np(phy.path)


#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################
geog.path <- here("inst", "extdata" , "geo_area_akodon.data")

# This is the example geography file for Hawaiian Psychotria
# (from Ree & Smith 2008)
geogfn = np(geog.path)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 4

numstates_from_numareas(numareas=5, 
                        maxareas=4, 
                        include_null_range=FALSE)


# run and save results of BioGeoBEARS -------------------------------------

  trfn = np(phy.path)
  geogfn = np(geog.path)
  
  #######################################################
  #######################################################
  # DEC AND DEC+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: The BioGeoBEARS "DEC" model is identical with 
  # the Lagrange DEC model, and should return identical
  # ML estimates of parameters, and the same 
  # log-likelihoods, for the same datasets.
  #
  # Ancestral state probabilities at nodes will be slightly 
  # different, since BioGeoBEARS is reporting the 
  # ancestral state probabilities under the global ML
  # model, and Lagrange is reporting ancestral state
  # probabilities after re-optimizing the likelihood
  # after fixing the state at each node. These will 
  # be similar, but not identical. See Matzke (2014),
  # Systematic Biology, for discussion.
  #
  # Also see Matzke (2014) for presentation of the 
  # DEC+J model.
  #######################################################
  #######################################################
  
  #######################################################
  #######################################################
  
  #######################################################
  # Run DEC
  #######################################################
  
  # Intitialize a default model (DEC model)
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  
  # Give BioGeoBEARS the location of the phylogeny Newick file
  BioGeoBEARS_run_object$trfn = trfn
  
  # Give BioGeoBEARS the location of the geography text file
  BioGeoBEARS_run_object$geogfn = geogfn
  
  # Input the maximum range size
  BioGeoBEARS_run_object$max_range_size = max_range_size
  
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  # 1. Here, un-comment ONLY the files you want to use.
  # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
  # 3. For example files see (a) extdata_dir, 
  #  or (b) http://phylo.wikidot.com/biogeobears#files
  #  and BioGeoBEARS Google Group posts for further hints)
  #
  # Uncomment files you wish to use in time-stratified analyses:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 2
  # (use more cores to speed it up; this requires
  # library(parallel) and/or library(snow). The package "parallel" 
  # is now default on Macs in R 3.0+, but apparently still 
  # has to be typed on some Windows machines. Note: apparently 
  # parallel works on Mac command-line R, but not R.app.
  # BioGeoBEARS checks for this and resets to 1
  # core with R.app)
  
  # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
  # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
  # but the results are imprecise and so I haven't explored it further.
  # In a Bayesian analysis, it might work OK, but the ML point estimates are
  # not identical.
  # Also, I have not implemented all functions to work with force_sparse=TRUE.
  # Volunteers are welcome to work on it!!
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  # BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC model
  # (nothing to do; defaults)
  
  # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
  BioGeoBEARS_run_object
  
  # This contains the model object
  BioGeoBEARS_run_object$BioGeoBEARS_model_object
  
  # This table contains the parameters of the model 
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
  
  # Run this to check inputs. Read the error messages if you get them!
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run DEC model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resDEC = res

  #######################################################
  # Run DEC+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 3
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  # BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run DECj model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resDECj = res
  
  
  #######################################################
  #######################################################
  # DIVALIKE AND DIVALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
  # Ronquist (1997)'s parsimony DIVA. It is a likelihood
  # interpretation of DIVA, constructed by modelling DIVA's
  # processes the way DEC does, but only allowing the 
  # processes DIVA allows (widespread vicariance: yes; subset
  # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
  #
  # DIVALIKE is a likelihood interpretation of parsimony
  # DIVA, and it is "like DIVA" -- similar to, but not
  # identical to, parsimony DIVA.
  #
  # I thus now call the model "DIVALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run DIVALIKE
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 3
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run DIVALIKE model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resDIVALIKE = res
  
  #######################################################
  # Run DIVALIKE+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 3
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  estart <- ifelse(estart < BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"], 
                   BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"],
                   estart)
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run DIVALIKEj model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resDIVALIKEj = res
  
  
  #######################################################
  #######################################################
  # BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
  # not identical with the full Bayesian model implemented 
  # in the "BayArea" program of Landis et al. (2013). 
  #
  # Instead, this is a simplified likelihood interpretation
  # of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
  # "d" and "e" work like they do in the DEC model of Lagrange 
  # (and BioGeoBEARS), and then BayArea's cladogenesis assumption
  # (which is that nothing in particular happens at cladogenesis) is 
  # replicated by BioGeoBEARS.
  #
  # This leaves out 3 important things that are in BayArea:
  # 1. Distance dependence (you can add this with a distances 
  #    matrix + the "x" parameter in BioGeoBEARS, however)
  # 2. A correction for disallowing "e" events that drive
  #    a species extinct (a null geographic range)
  # 3. The neat Bayesian sampling of histories, which allows
  #    analyses on large numbers of areas.
  #
  # The main purpose of having a "BAYAREALIKE" model is 
  # to test the importance of the cladogenesis model on 
  # particular datasets. Does it help or hurt the data 
  # likelihood if there is no special cladogenesis process?
  # 
  # BAYAREALIKE is a likelihood interpretation of BayArea,
  # and it is "like BayArea" -- similar to, but not
  # identical to, Bayesian BayArea.
  # I thus now call the model "BAYAREALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run BAYAREALIKE
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 3
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  # BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE model
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # Check the inputs
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run BAYAREALIKE model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resBAYAREALIKE = res
  
  
  #######################################################
  # Run BAYAREALIKE+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "data/BioGeoBEARS/timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "data/BioGeoBEARS/areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = "GenSA"
  BioGeoBEARS_run_object$num_cores_to_use = 3
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  # BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run BAYAREALIKE model and save results
  res = bears_optim_run(BioGeoBEARS_run_object)
  resBAYAREALIKEj = res  

  #########################################################################
  #########################################################################
  #########################################################################
  #########################################################################
  # 
  # CALCULATE SUMMARY STATISTICS TO COMPARE
  # DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  # 
  #########################################################################
  #########################################################################
  #########################################################################
  #########################################################################
  
  #########################################################################
  #########################################################################
  # REQUIRED READING:
  #
  # Practical advice / notes / basic principles on statistical model 
  #    comparison in general, and in BioGeoBEARS:
  # http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
  #########################################################################
  #########################################################################
  
  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL
  
  #######################################################
  # Statistics -- DEC vs. DEC+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DEC, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DEC+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
  # confer the same likelihood on the data. See: Brian O'Meara's webpage:
  # http://www.brianomeara.info/tutorials/aic
  # ...for an intro to LRT, AIC, and AICc
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- DIVALIKE vs. DIVALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DIVALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #########################################################################
  # ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  #########################################################################
  teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
  teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
  row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
  restable = put_jcol_after_ecol(restable)
  
  ########################################################
  ## Save the results tables for later -- check for e.g.
  ## convergence issues
  ########################################################
  #
  ## Loads to "restable"
  #save(restable, file="restable_v1.Rdata")
  #load(file="restable_v1.Rdata")
  #
  ## Loads to "teststable"
  #save(teststable, file="teststable_v1.Rdata")
  #load(file="teststable_v1.Rdata")
  #
  ## Also save to text files
  #write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
  #write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")
  
  #######################################################
  # Model weights of all six models
  #######################################################
  restable2 = restable
  
  # With AICs:
  AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable = cbind(restable, AICtable)
  restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
  restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
  restable_AIC_rellike
  
  # With AICcs -- factors in sample size
  tr <- read.tree(trfn)
  samplesize = length(tr$tip.label)
  AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
  restable2 = cbind(restable2, AICtable)
  restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
  restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
  restable_AICc_rellike
  
  
  list.results.BioGeoBears <- list(resDEC = resDEC, 
                                   resDECj = resDECj,
                                   resDIVALIKE  = resDIVALIKE ,
                                   resDIVALIKEj = resDIVALIKEj,
                                   resBAYAREALIKE  = resBAYAREALIKE ,
                                   resBAYAREALIKEj = resBAYAREALIKEj,
                                   table_AIC = restable_AIC_rellike ,
                                   table_AICc = restable_AICc_rellike)  

  list.results.BioGeoBears$table_AICc %>% 
    dplyr::arrange(AICc) %>% 
    mutate(delta_AICc =  AICc - min(AICc))
  
  save(resDEC, file =  here("inst", "extdata", "resDEC_akodon.RData"))
  
  cat("The best models are 'resDEC', 'resBAYAREALIKEj', and 'resBAYAREALIKE'")
  