#!/usr/bin/bash
#SBATCH -J squire_prepare                  # A single job name for the array
#SBATCH -n 10                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 100000                # Memory request (16Gb)
#SBATCH -t 10-00:00                  # Maximum execution time (D-HH:MM)
#SBATCH -p bigmem

echo 'squire Fetch at:' `date`
squire Fetch --build hg38 --fasta --rmsk  --chrom_info  --index  --gene --pthreads 10
echo 'Complete Fetch at:' `date`
#*******************

