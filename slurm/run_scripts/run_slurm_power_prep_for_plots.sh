#!/bin/bash

cp ~/ACEBounds/scripts/power/power_prep_for_plots.R /workspace/rtrane/ACEBounds/scripts/power/.

sbatch ~/ACEBounds/slurm/slurm_scripts/slurm_power_prep_for_plots.sh
