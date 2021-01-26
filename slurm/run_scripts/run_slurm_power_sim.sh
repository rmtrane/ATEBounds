#!/bin/bash

cp power_sims_only.R /workspace/rtrane/ACEBounds/.

sbatch slurm_power_sim.sh
