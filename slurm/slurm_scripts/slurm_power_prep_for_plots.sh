#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rtrane@wisc.edu
#SBATCH -o prep_for_plots.out
#SBATCH -e prep_for_plots.error
#SBATCH -D /workspace/rtrane/ACEBounds
#SBATCH -J prep_for_plots
#SBATCH -t 72:00:00
#SBATCH -p long
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --array=56 # 1-84

# Make sure R knows where I keep my packages
export R_LIBS=/workspace/rtrane/Rpackages

module load R/R-4.0.1

mkdir V8_$SLURM_ARRAY_TASK_ID

R CMD BATCH --no-save --no-restore "--args $SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID" power_prep_for_plots.R Routs/power_prep_for_plots_$SLURM_ARRAY_TASK_ID.Rout

rm -rf V8_$SLURM_ARRAY_TASK_ID
