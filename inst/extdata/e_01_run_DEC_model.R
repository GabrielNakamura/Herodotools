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

phy.path <-  here("examples", "akodon.new")

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
geog.path <- here("examples", "geo_area_akodon.data")

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
  
  

