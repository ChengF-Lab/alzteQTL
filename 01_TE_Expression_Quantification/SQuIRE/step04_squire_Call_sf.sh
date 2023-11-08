#!/usr/bin/bash
#SBATCH -J squire_Call              # A single job name for the array
#SBATCH -n 10                       # Number of cores
#SBATCH --mem 80000                 # Memory request (16Gb)
#SBATCH -t 10-00:00                 # Maximum execution time (D-HH:MM)
#SBATCH -p defq

treat=`cat $1`
control=`cat $2`
dir=`pwd`

echo "squire start call step at:" `date`

squire Call -1 $treat -2 $control -A AD -B HC -o $dir/squire_Call_sf -p 10 -s True -N sf -f pdf

echo "squire end call step at:" `date`


