# Data analysis/figures reconstruction from O'Brien et al. 2023

## Overview
This repository contains scripts used to produce the results and figures in [O'Brien et al. 2023](www.google.com).
To run the scripts, you will need to download the data from https://doi.org/10.48610/f3850b0 and follow the steps below.

## Requirements
You will need at least 8GB of RAM and a Linux operating system to run the analysis.
You will also need to build and install our [ODELandscaper software](https://github.com/nobrien97/odeLandscape/tree/main/ODESolver) to calculate fitness effects in NAR models.
The analysis was done in R 4.3.1, using the following packages (with versions):
- dplyr (1.1.2)
- ggplot2 (3.4.2)
- tibble (3.2.1)
- tidyr (1.3.0)
- cowplot (1.1.1)
- ggridges (0.5.4)
- ggpmisc (0.5.3)
- deSolve (1.35)
- DescTools (0.99.49)
- paletteer (1.5.0)
- latex2exp (0.9.6)
- ggraph (2.1.0)
- igraph (1.4.2)
- fitdistrplus (1.1.8)
- GenSA (1.1.8)
- future (1.32.0)
- doParallel (1.0.17)
- foreach (1.5.2)

## Step 1
Clone this repo (or download the release) and the dataset. Extract the dataset to a convenient directory.
In `runAnalysis.R`, replace the paths on lines 8 and 9 with the paths to where you saved this repo, and to
where you saved the dataset.

## Step 2
Run runAnalysis.R 
This will save the figures from the paper into the main repo directory (`NARAdaptiveWalk2023/`).
Stats will be printed to the R console -- look at the source in stats.R for more information.


## Running SLiM simulations
You can find the SLiM scripts used in the paper in the `./SLiM/` folder, along with a list of random number seeds used. We used a custom version of SLiM available at [](https://github.com/nobrien97/SLiM/releases/tag/AdaptiveWalks2023). Note that to run the scripts, you will have to change a few variables -- the `wd` constant on line 12, and `moveDir` (directory to move all data to after the simulation is done) on line 80 (`indTrack_add.slim`) and line 84 (`indTrack_net.slim`).

To run the simulation, it is easiest to use the command line. An example PBS-style script is below:
```
#!/bin/bash -l
#PBS -l walltime=24:00:00
#PBS -l ncpus=1440
#PBS -l mem=2000GB
  
cd /PATH/TO/NARAdaptiveWalk2023/SLiM
SAVEDIR=/PATH/TO/SAVED/FILES

# Run the SLiM scripts
for seed in seeds.csv; do
  echo "Running with seed == " $(seed):
  slim -s $(seed) indTrack_add.slim &
  slim -s $(seed) indTrack_net.slim &
  echo 
done

echo "All jobs finished, moving output..."

# Combine output into a single file
cd /PATH/TO/moveDir # set up in .slim files
SAVEDIR=/PATH/TO/SAVE/COMBINED/FILES

cat ./slim_muts* >> $SAVEDIR/slim_muts.csv
cat ./slim_qg* >> $SAVEDIR/slim_qg.csv
cat ./slim_locusHo* >> $SAVEDIR/slim_locusHo.csv

# Saves population states after burn-in if you want to run them again
mkdir -p $SAVEDIR/popstates
mv ./slim_popstate* $SAVEDIR/popstates

# Delete loose files
find -regex ".*[0-9]*_*[0-9].csv+" -delete

```
A smarter way to run these scripts is to use a high-throughput approach like MPIrun or job arrays to maximise
CPU usage, but the implementation depends on the HPC and its available software.

Note that each NAR simulation takes around 13 hours to run on an Intel Xeon Platinum 8274 3.2GHz core. Additive simulations are slightly faster, at around 9 hours. We ran the simulations on the NCI Gadi HPC system, using a 
high-throughput design to massively parallelise the 5760 simulations.