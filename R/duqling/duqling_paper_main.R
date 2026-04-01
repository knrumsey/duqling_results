library(duqling)
library(tidyverse)

# Choose file extension for figures
extension <- "png"

# Warnings are normal from these scripts.

# Section 4 Figures
source("R/duqling/make_figs_main.R")
#source("R/make_figs_SM.R")

# Section 5 Figures
source("R/duqling/no_free_lunch.R")
source("R/duqling/tuning_smallnug.R")
source("R/duqling/clustering.R")
source("R/duqling/clustering_data.R")
