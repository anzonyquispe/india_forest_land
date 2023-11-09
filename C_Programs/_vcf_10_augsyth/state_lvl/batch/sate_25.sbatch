#!/bin/bash
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=90G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-user=futurolos9@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=sate_25.log

module purge
R CMD BATCH /scratch/gpfs/ar8787/groupdata2/india_forest_land/C_Programs/_vcf_10_augsyth/state_lvl/batch/sate_25.R
