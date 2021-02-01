#!/bin/bash

mkdir -p ../data/power
mkdir -p ../power_sims_out/

for i in {1..840}; do
    R CMD BATCH --no-save --no-restore "--args $i" power/power_sims_only.R ../power_sims_out/power_sims_$i.Rout
done

for i in {1..84}; do
    R CMD BATCH --no-save --no-restore "--args 1 $i" power/power_prep_for_plots.R ../power_sims_out/power_prep_for_plots_$i.Rout
done

R CMD BATCH --no-save --no-restore power/power_plots_and_results.R ../power_sims_out/power_plots_an_results.R
