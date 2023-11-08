#!/usr/bin/bash
#SBATCH -J squire_Map               # A single job name for the array
#SBATCH -n 10                       # Number of cores
#SBATCH --mem 80000                 # Memory request (16Gb)
#SBATCH -t 10-00:00                 # Maximum execution time (D-HH:MM)
#SBATCH -p defq

raw_dir=$1
sample=$2

dir=`pwd`
cd $dir

cd 01_Clean_data   #### You should first mkdir 01_Clean_data in your output dir  ####

echo "$sample quality control using fastp start at:" `date`
   fastp -i $raw_dir/$sample.fq.gz -o $sample.fq.gz -z 4 -q 20 -u 40 -n 5 -j $sample.json -h $sample.html
echo "$sample quality control using fastp end at:" `date`
cd ../

echo "$sample start squire Map at:" `date`
   squire Map -1 01_Clean_data/$sample.fq.gz -r 50 -p 15 -n $sample -f squire_fetch
echo "$sample end squire Map at:" `date`

echo "$sample start squire Count at:" `date`
   squire Count -r 50 -p 15 -n $sample -c squire_clean -f squire_fetch 
echo "$sample end squire Count at:" `date`


