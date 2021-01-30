#!/bin/bash

cp ~/ACEBounds/scripts/power/* /workspace/rtrane/ACEBounds/scripts/power/.

sbatch ~/ACEBounds/slurm/slurm_scripts/slurm_power_sim.sh
