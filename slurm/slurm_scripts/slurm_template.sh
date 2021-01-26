#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rtrane@wisc.edu
#SBATCH -o small_from_sparse.out
#SBATCH -e small_from_sparse.error
#SBATCH -D /workspace/rtrane/BMI826/
#SBATCH -J small_from_sparse
#SBATCH -t 72:00:00
#SBATCH -p long
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=200M

module load python/python2.7.14

cd BETS/code

bash run_BETS.sh small_from_large_from_sparse
