# Run O'Brien et al. 2023 analysis: 
#   format data, run mutation experiment, statistics, figures
# Author: Nick O'Brien
# Last updated: 2023/08/08

# Replace these paths with the path to where you saved the 
# NarAdaptiveWalk2023 repo and where you saved the dataset
repoPath <- "/path/to/NARAdaptiveWalk2023"
dataPath <- "/path/to/data/"
setwd(repoPath)

# Load functions
source("./R/helperFunctionsAndSetup.R")

# Setup data
source("./R/wrangle_data.R")

# Run mutation screen experiment
source("./R/mutationScreenExp.R")

# Load plotting functions
source("./R/figureFunctionsAndSetup.R")

# Create figures
source("./R/figures.R")

# Run stats
source("./R/stats.R")
