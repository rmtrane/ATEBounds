#!/bin/bash

R CMD BATCH --no-save --no-restore multiple_IVs/multiple_IVs_sims.R multiple_IVs/multiple_IVs_sims.R 
R CMD BATCH --no-save --no-restore multiple_IVs/multiple_IVs_prep.R multiple_IVs/multiple_IVs_prep.R
R CMD BATCH --no-save --no-restore multiple_IVs/multiple_IVs_results_and_figures.R multiple_IVs/multiple_IVs_results_and_figures.R 
