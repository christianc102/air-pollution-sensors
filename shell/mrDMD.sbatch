#!/bin/bash
#SBATCH -J mrDMD
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-10:00
#SBATCH -p serial_requeue
#SBATCH --mem=100000
#SBATCH --mail-type=END

module load matlab/R2021a-fasrc01
matlab -nodesktop -nosplash -nodisplay < "mrDMD_second.m"

exit 0
