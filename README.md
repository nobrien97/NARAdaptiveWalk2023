# Data analysis/figures reconstruction from O'Brien et al. 2023

## Overview
This repository contains scripts used to produce the results and figures in [O'Brien et al. 2023](www.google.com).
To run the scripts, you will need to download the data from [Dryad](www.google.com) and follow the steps below.

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



