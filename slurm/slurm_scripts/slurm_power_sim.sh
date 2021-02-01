#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rtrane@wisc.edu
#SBATCH -o power_sims.out
#SBATCH -e power_sims.error
#SBATCH -D /workspace/rtrane/ACEBounds
#SBATCH -J power_sims
#SBATCH -t 72:00:00
#SBATCH -p long
#SBATCH --mem-per-cpu=8000M
#SBATCH --array=1-840

# Make sure R knows where I keep my packages
export R_LIBS=/workspace/rtrane/Rpackages

module load R/R-4.0.1

mkdir V8_$SLURM_ARRAY_TASK_ID

R CMD BATCH --no-save --no-restore "--args $SLURM_ARRAY_TASK_ID" scripts/power/power_sims_only.R power_sims_out/power_sims_$SLURM_ARRAY_TASK_ID.Rout

rm -rf V8_$SLURM_ARRAY_TASK_ID
