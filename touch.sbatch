#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-8:00
#SBATCH -p test
#SBATCH --mem=8000
#SBATCH --mail-type=END

find . -type f -exec touch {} +

# Exit normally
exit 0
#EOC
