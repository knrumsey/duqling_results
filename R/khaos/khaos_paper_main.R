library(duqling)
library(tidyverse)
library(khaos)

rerun_simulations <- FALSE # We can just load these instead
rerun_sobol <- FALSE       # Do you want to recreate Sobol figures?

#==================================
#       RUN SIMULATION STUDY
#==================================
source("R/khaos/khaos_sim_study.R")

# Or just load in the data:
load("data/results_khaos.Rda")

#==================================
#       MAKE FIGURES
#==================================
source("R/khaos/khaos_figures.R")

# Takes slightly longer to run these:
if(rerun_sobol){
  source("R/khaos/khaos_sobol.R")
}


