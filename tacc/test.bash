#!/bin/bash
#SBATCH -J future-test          # Job name
#SBATCH -o future.o%j           # Name of stdout output file (%j expands to jobId)
#SBATCH -e future.o%j           # Name of stderr output file(%j expands to jobId)
#SBATCH -p skx-dev			    # Submit to the 'normal' or 'development' queue
#SBATCH -N 1                    # Total number of nodes
#SBATCH -n 48                  	# Total number of mpi tasks requested
#SBATCH -t 2:00:00             	# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load custom R module
eval "$(conda shell.bash hook)"
conda activate custom-R

Rscript TACC-script-multicore.R