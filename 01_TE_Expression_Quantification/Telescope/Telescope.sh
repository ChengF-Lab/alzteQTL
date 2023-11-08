#!/usr/bin/bash

input="squire_map"   ##### your bam file folder
sample=`cat $1`      ##### your sample list

for sam in $sample
	do
		echo "$sam START at date: "`date`
		
		telescope assign $input/$sam\.bam total_TE.gtf --exp_tag $sam
		telescope resume ./$sam-checkpoint.npz --exp_tag $sam

		echo "$sam END at date: "`date`
done


